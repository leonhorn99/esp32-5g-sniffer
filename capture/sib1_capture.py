"""Capture full bursts timed to a cell's SIB1 slots, for matlab/sib1_offline.m.

  .venv_local/bin/python -m capture.sib1_capture --port /dev/ttyACM0 --freq 2670 --ssb 2667.75 --pci 317 --sfn-mod 16 --slot 3
  .venv_local/bin/python -m capture.sib1_capture --port /dev/ttyACM0 --freq 2130 --ssb 2125.45 --pci 165
  .venv_local/bin/python -m capture.sib1_capture --port /dev/ttyACM0 --freq 2670 --ssb 2667.75 --pci 317 --sfn-mod 64 --sfn-set 0 2 4 6 --label si
  .venv_local/bin/python -m capture.sib1_capture --port /dev/ttyACM0 --freq 2130 --ssb 2125.45 --pci 165 --sfn-mod 8 --sfn-set 0 --slots 0 1 2 3 4 5 6 7 8 9 --label si
  then: matlab -batch "addpath('matlab'); sib1_offline('data/<session>', 'FreqMHz', 2670, 'SSBMHz', 2667.75, 'PCI', 317)"

A burst (1.64 ms at 20 MS/s) that holds the SSB never holds the slots of SIB1's control channel, and
untimed bursts hit such a slot only by chance (5 of 60 bursts on 23 Nov; with SIB1 only every 160 ms,
2 of 60). So:
1. SSB: bursts timed to cover the 20 ms SSB period (cell_scan.Coverage) until one holds the cell's
   SSB with a decoded MIB (NumPy receiver at the SSB's frequency, the PCI's N_ID2 only).
2. Frame clock: board time of frame k = a + b k, b from the CFO (the board's crystal: 10 ms x
   (1 - CFO / f_SSB)), a from the SSB anchors (refitted with every new one).
3. From the MIB (capture/nr_mib.py, TS 38.213 Tables 13-1, 13-11): the Type0-PDCCH slots n0, n0+1
   of even or odd frames, with --sfn-mod M only frames whose SFN is a multiple of M.
4. Every further burst is delayed so that the next such slot lies in its middle; every --anchor-every-th
   burst is put on an SSB instead (fresh anchors for the clock and for sib1_offline's own fit). The
   bursts' spacing on the board's clock is measured as it goes; it jitters by about +-1 ms against a
   0.64 ms window, so not every burst holds a whole slot. Until --hits bursts do.
All bursts are recorded as a session in data/. Receive only.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from .cell_scan import Coverage, cfo_window
from .nr_mib import sib1_occasions

ROOT = Path(__file__).resolve().parents[1]
FRAME_US = 10000.0
SLOT_US = 1000.0                                  # 15 kHz
SSB_US = 4 * 1e6 / 15e3 * (1 + 144 / 2048) + 2   # 4 OFDM symbols with CP (286 us) + margin
MARGIN_US = 20.0                                  # a slot counts as held if it lies this far inside the burst


def ssb_anchor(x, fs: float, f_lo: float, f_ssb: float, pci: int, cfo_c: float, cfo_span: float,
               t0_us: float, scale: float = 1.0) -> dict | None:
    """The cell's SSB in one burst, if its MIB decodes: {t_frame_us, sfn, cfo, issb, mib}."""
    from . import nr_ssb as nr
    cfos = round(cfo_c / 250) * 250 + np.arange(-int(cfo_span // 5000), int(cfo_span // 5000) + 1) * 5000
    for p in nr.pss_search_raster(x, fs, [f_ssb - f_lo], cfos, nid2s=(pci % 3,)):
        r = nr.receive(x, fs, p["t_pss"], f_ssb, f_lo, cfo=p["cfo"], nid2=p["nid2"], t_first_us=t0_us, scale=scale)
        if r.get("ok") and r.get("pci") == pci and r.get("t_frame_us") is not None:
            return {"t_frame_us": r["t_frame_us"], "sfn": r["mib"]["sfn"], "cfo": r["cfo"], "issb": r["issb"],
                    "mib": r["mib"]}
    return None


class FrameClock:
    """Board time (us) of the start of frame k (SFN unwrapped): a + b k."""

    def __init__(self, us_per_frame: float):
        self.b = us_per_frame
        self.anchors: list[tuple[int, float]] = []
        self.a = None

    def add(self, t_frame_us: float, sfn: int):
        if self.anchors:
            k0, t0 = self.anchors[0]
            approx = k0 + (t_frame_us - t0) / self.b
            k = sfn + 1024 * round((approx - sfn) / 1024)
        else:
            k = sfn
        self.anchors.append((k, t_frame_us))
        ks = np.array([k for k, _ in self.anchors], float)
        ts = np.array([t for _, t in self.anchors])
        self.a = float(np.mean(ts - self.b * ks))      # rate fixed by the CFO; offset from all anchors
        return k

    def t(self, k: float) -> float:
        return self.a + self.b * k

    def frame_at(self, t_us: float) -> float:
        return (t_us - self.a) / self.b


def sib1_targets(clock: FrameClock, mib: dict, issb: int, after_us: float, dur_us: float, sfn_mod: int = 1,
                 n: int = 4, only_slot: int | None = None, sfn_set=None, slots=None) -> list[tuple[float, int, int]]:
    """The next n burst start times (board us) that centre a Type0-PDCCH slot (n0 or n0+1, or only_slot)
    in a burst of dur_us, with (SFN, slot), from after_us on. Frames: SFN % sfn_mod == 0, or with
    sfn_set SFN % sfn_mod in sfn_set (e.g. an SI window: mod 64 in 0, 2, 4, 6)."""
    if slots is None:                          # the Type0-PDCCH slots of SIB1 (n0, n0 + 1 of even/odd frames)
        occ = sib1_occasions(mib, issb)
        if occ is None:
            raise ValueError("the MIB announces no SIB1 control channel this tool can find (k_SSB >= 24, SCS or table)")
        parity, slot_list = occ["frame_parity"], occ["slots"]
    else:                                      # explicit slots of any frame (e.g. a search space for other SI)
        parity, slot_list = None, list(slots)
    out, k = [], int(np.floor(clock.frame_at(after_us))) - 1
    while len(out) < n:
        k += 1
        sfn = k % 1024
        if (parity is not None and sfn % 2 != parity) or \
                (sfn % sfn_mod not in sfn_set if sfn_set else sfn % sfn_mod):
            continue
        for s in slot_list:
            if only_slot is not None and s != only_slot:
                continue
            start = clock.t(k + s / 10) - (dur_us - SLOT_US) / 2
            if start >= after_us:
                out.append((start, sfn, s))
    return out[:n]


def ssb_target(clock: FrameClock, t_ssb_in_frame_us: float, after_us: float, dur_us: float) -> float:
    """The next burst start (board us) that centres the cell's SSB (period 20 ms) in a burst."""
    k = int(np.floor(clock.frame_at(after_us))) - 2
    while True:
        k += 1
        start = clock.t(k) + t_ssb_in_frame_us - (dur_us - SSB_US) / 2
        if start >= after_us and k % 2 == clock.anchors[0][0] % 2:     # SSB every 2nd frame
            return start


def holds_slot(clock: FrameClock, t0_us: float, dur_us: float, mib: dict, issb: int, sfn_mod: int = 1,
               only_slot: int | None = None, sfn_set=None, slots=None):
    """(SFN, slot) of a target slot that lies whole in the burst, or None."""
    for start, sfn, s in sib1_targets(clock, mib, issb, t0_us - dur_us, dur_us, sfn_mod, n=12, only_slot=only_slot,
                                      sfn_set=sfn_set, slots=slots):
        ts = start + (dur_us - SLOT_US) / 2
        if t0_us + MARGIN_US <= ts and ts + SLOT_US <= t0_us + dur_us - MARGIN_US:
            return sfn, s
    return None


MAX_OFFSET_MHZ = 9.5   # CORESET#0 and SSB at most this far from the LO (20 MS/s: +-10 MHz; the C5 passes +-8 well)
DC_GAP_MHZ = 0.1       # the LO (the receiver's DC) at least this far outside CORESET#0


def lo_for_coreset(mib: dict, f_ssb_hz: float, lo_mhz: int | None = None) -> int:
    """An LO (whole MHz) whose 20 MS/s burst holds the cell's CORESET#0 (and so SIB1's PDSCH, which lies in
    the same RBs) and its SSB, with the DC outside CORESET#0: the given LO if it does, else the tunable
    one with the smallest largest offset. PCI 165 (CORESET#0 2120.68-2129.32 MHz): 2130; PCI 317: 2670."""
    from .nr_mib import coreset0_band_hz
    from .ssb_sources import lo_tunable
    band = coreset0_band_hz(mib, f_ssb_hz)
    if band is None:
        raise ValueError("no CORESET#0 in this MIB (k_SSB >= 24 or a table this tool does not have)")
    lo_b, hi_b = band[0] / 1e6, band[1] / 1e6
    edges = (lo_b, hi_b, f_ssb_hz / 1e6 - 1.8, f_ssb_hz / 1e6 + 1.8)

    def worst(lo):
        if lo_b - DC_GAP_MHZ < lo < hi_b + DC_GAP_MHZ or not lo_tunable(lo):
            return None
        w = max(abs(e - lo) for e in edges)
        return w if w <= MAX_OFFSET_MHZ else None

    if lo_mhz is not None and worst(int(lo_mhz)) is not None:
        return int(lo_mhz)
    cands = [(worst(lo), lo) for lo in range(int(lo_b) - 10, int(hi_b) + 12)]
    cands = [c for c in cands if c[0] is not None]
    if not cands:
        raise ValueError(f"no tunable LO holds CORESET#0 {lo_b:.2f}-{hi_b:.2f} MHz with the DC outside it")
    return min(cands)[1]


TICK_US = 1000.0       # a capture starts on the next USB full-speed frame (1 ms) after its command arrives
LEAD_US = 3000.0       # least time from now to a planned start (host sleep, USB, command parsing)


class TickModel:
    """Where captures can start: the 1 ms USB ticks, in board time, and the offset D between the host's
    clock (when a command is sent) and the board's: start = first tick >= host time + D. Measured on
    26 Nov: the tick grid drifts against the cell's slots by about 5 us/s, so for minutes at a time no
    tick lets a 1.64 ms burst hold a whole 1 ms slot, and then for minutes every slot can be hit."""

    def __init__(self):
        self.pairs: list[tuple[float, float]] = []          # (host us at the command, board start us)

    def add(self, host_us: float, t0_us: float):
        self.pairs.append((host_us, t0_us))
        del self.pairs[:-20]

    def ready(self) -> bool:
        return len(self.pairs) >= 3

    def offset(self) -> float:
        lo = max(t0 - TICK_US - h for h, t0 in self.pairs)
        hi = min(t0 - h for h, t0 in self.pairs)
        return (lo + hi) / 2 if lo <= hi else float(np.median([t0 - h - TICK_US / 2 for h, t0 in self.pairs]))

    def tick_in(self, lo_us: float, hi_us: float) -> float | None:
        """The first tick in [lo_us, hi_us] (board us), or None."""
        last = self.pairs[-1][1]
        t = last + np.ceil((lo_us - last) / TICK_US) * TICK_US
        return float(t) if t <= hi_us else None


def run(sdr, a, rate_code: int, writer=None, say=print, sleep=time.sleep, now=time.monotonic) -> dict:
    """The capture loop (board already tuned). Returns counts, anchors and the slots hit."""
    from .ssb_sources import iq_scale
    f_lo, f_ssb = a.freq * 1e6, a.ssb * 1e6
    cfo_c, cfo_span = cfo_window(a.ssb, a.ppm, 30e3)
    cov, clock, mib, issb, t_ssb = Coverage(), None, None, None, None
    last = None; gaps = []; delay = 0.0; dur = None
    hits, n, n_anchor_tries, per_class = [], 0, 0, {}
    want_ssb, waits = False, 0
    ticks = TickModel()
    while n < a.max_bursts and len(hits) < a.hits:
        if clock is not None and dur and ticks.ready():
            # plan on the tick grid: the tick that puts the target whole into the burst
            D = ticks.offset()
            earliest = now() * 1e6 + D + LEAD_US
            start = None
            if not want_ssb:
                cands = sib1_targets(clock, mib, issb, earliest, dur, a.sfn_mod, n=24, only_slot=a.slot,
                                     sfn_set=getattr(a, "sfn_set", None), slots=getattr(a, "slots", None))
                cands.sort(key=lambda c: (per_class.get((c[1] % a.sfn_mod, c[2]), 0), c[0]))
                for c0, sfn, sl in cands:
                    slot_start = c0 + (dur - SLOT_US) / 2
                    t = ticks.tick_in(max(earliest, slot_start + SLOT_US + MARGIN_US - dur), slot_start - MARGIN_US)
                    if t is not None:
                        start = t
                        break
                if start is None:                       # no tick fits a slot now: an SSB anchor, then wait
                    want_ssb = True
                    waits += 1
                    if waits % 20 == 1:
                        say("  the USB tick grid misses the SIB1 slots for now (it drifts about 5 us/s against "
                            "them); taking SSB anchors meanwhile")
                    sleep(0.5)
                    D = ticks.offset()
                    earliest = now() * 1e6 + D + LEAD_US
            if want_ssb:
                k0 = int(np.floor(clock.frame_at(earliest)))
                for k in range(k0, k0 + 6):
                    if k % 2 != clock.anchors[0][0] % 2:
                        continue
                    s0 = clock.t(k) + t_ssb
                    start = ticks.tick_in(max(earliest, s0 + SSB_US + MARGIN_US - dur), s0 - MARGIN_US)
                    if start is not None:
                        break
            if start is not None:
                sleep(max(0.0, (start - D - TICK_US / 2) / 1e6 - now()))
        elif last is not None and gaps and dur:
            # before the clock (or the tick model) exists: bursts timed to cover the SSB period
            spacing = float(np.median(gaps[-8:]))
            pred = last + spacing
            delay = (cov.next_start(1, dur) - pred) % cov.period if clock is None else 0.0
            sleep(delay / 1e6)
        else:
            delay = 0.0
        h_cmd = now() * 1e6
        b = sdr.capture(sdr.max_samples, rate_code, bits=a.bits)
        if last is not None:
            gaps.append(b.t0_us - last - delay)
        last = b.t0_us
        ticks.add(h_cmd, b.t0_us)
        if not b.crc_ok:
            continue
        n += 1
        if writer:
            writer.add(b, freq_mhz=a.freq, gain=a.gain, bandwidth_mhz=None)
        x = b.iq()
        dur = len(x) / b.fs * 1e6
        if clock is None or want_ssb:
            n_anchor_tries += 1
            anc = ssb_anchor(x, b.fs, f_lo, f_ssb, a.pci, cfo_c, cfo_span, b.t0_us, iq_scale(b.fmt))
            if anc:
                if clock is None:
                    ppm = -anc["cfo"] / a.ssb                    # Hz / MHz
                    clock = FrameClock(FRAME_US * (1 + ppm * 1e-6))
                    cfo_c, cfo_span = cfo_window(a.ssb, ppm, 30e3)
                    mib, issb = anc["mib"], anc["issb"]
                    t_ssb = _ssb_offset_us(issb, mib["hrf"])     # the SSB's start within its frame
                    say(f"SSB of PCI {a.pci} in burst {n}: SFN {anc['sfn']}, SSB index {issb}, CFO "
                        f"{anc['cfo'] / 1e3:+.2f} kHz ({ppm:+.2f} ppm), pdcch-ConfigSIB1 {mib['pdcch_config_sib1']}")
                clock.add(anc["t_frame_us"], anc["sfn"])
            elif clock is None:
                cov.add(b.t0_us, dur)
            want_ssb = False
        else:
            h = holds_slot(clock, b.t0_us, dur, mib, issb, a.sfn_mod, a.slot, getattr(a, "sfn_set", None),
                           getattr(a, "slots", None))
            if h:
                hits.append((n, *h))
                per_class[(h[0] % a.sfn_mod, h[1])] = per_class.get((h[0] % a.sfn_mod, h[1]), 0) + 1
                say(f"  burst {n}: slot {h[1]} of SFN {h[0]} inside ({len(hits)}/{a.hits})")
            want_ssb = (n % a.anchor_every == 0)
    return {"bursts": n, "hits": hits, "anchors": len(clock.anchors) if clock else 0, "mib": mib, "issb": issb,
            "anchor_tries": n_anchor_tries, "waits": waits}


def _ssb_offset_us(issb: int, hrf: int) -> float:
    """Start of SSB i (15 kHz case A) within its frame: first symbols 2, 8 + 14 n (TS 38.213 4.1)."""
    first = [2, 8, 16, 22][issb]
    slot, sym = divmod(first, 14)
    t_sym = 1e6 / 15e3 * (1 + 144 / 2048)
    return 5000 * hrf + slot * SLOT_US + 16 / 2048 * 1e6 / 15e3 + sym * t_sym


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--port", required=True)
    p.add_argument("--freq", type=int, required=True, help="LO in whole MHz; the CORESET#0 must fit into +-8 MHz")
    p.add_argument("--ssb", type=float, required=True, help="SSB centre (SS_REF) in MHz")
    p.add_argument("--pci", type=int, required=True)
    p.add_argument("--hits", type=int, default=20, help="bursts that hold a whole SIB1 slot")
    p.add_argument("--sfn-mod", type=int, default=1, help="only frames whose SFN is a multiple of this (16: SIB1 "
                                                          "once per 160 ms, as PCI 317 on 23 Nov)")
    p.add_argument("--sfn-set", type=int, nargs="*", default=None,
                   help="frames with SFN %% --sfn-mod in this set, e.g. --sfn-mod 64 --sfn-set 0 2 4 6 for the SI "
                        "window of SIB2-5 of PCI 317 (SI messages use the SIB1 occasions when "
                        "searchSpaceOtherSystemInformation is 0)")
    p.add_argument("--slots", type=int, nargs="+", default=None,
                   help="target these slots of every selected frame instead of SIB1's Type0 slots, e.g. 0 1 ... 9 for "
                        "an SI window when searchSpaceOtherSystemInformation is not 0 (PCI 165: search space 1, "
                        "every slot; SIB2-4 in the frame with SFN %% 8 = 0)")
    p.add_argument("--slot", type=int, default=None, help="only this Type0-PDCCH slot (default both, n0 and n0+1; "
                                                          "PCI 317 used n0 = 3 on 23 Nov)")
    p.add_argument("--anchor-every", type=int, default=8, help="every n-th burst on an SSB (frame clock)")
    p.add_argument("--max-bursts", type=int, default=600)
    p.add_argument("--ppm", type=float, default=None, help="board crystal, if known (narrows the first SSB search)")
    p.add_argument("--rate", type=float, default=20, help="MS/s")
    p.add_argument("--bits", type=int, choices=(8, 10), default=8)
    p.add_argument("--gain", default="auto")
    p.add_argument("--label", default="")
    a = p.parse_args(argv)
    from .capture_bursts import _setup, pick_gain
    from .esp_sdr import EspSdr, rate_code_for
    from .session import SessionWriter
    sdr = EspSdr(a.port, timeout=3.0)
    device = sdr.connect()
    rate = rate_code_for(a.rate * 1e6)
    a.bandwidth, a.dccal = None, False
    if a.gain == "auto":
        a.gain, _ = pick_gain(sdr, a.freq, rate, 10, a)
        print(f"auto gain: {a.gain}")
    a.gain = int(a.gain)
    _setup(sdr, a.freq, a, a.gain)
    name = f"{time.strftime('%Y%m%d_%H%M%S')}_sib1_pci{a.pci}" + (f"_{a.label}" if a.label else "")
    writer = SessionWriter(ROOT / "data" / name, device,
                           {"tool": "sib1_capture", "freq_mhz": [a.freq], "ssb_mhz": a.ssb, "pci": a.pci,
                            "rate_msps": a.rate, "bits": a.bits, "gain": a.gain, "sfn_mod": a.sfn_mod})
    t0 = time.monotonic()
    try:
        res = run(sdr, a, rate, writer)
    finally:
        writer.close()
        sdr.close()
    print(f"{res['bursts']} bursts in {time.monotonic() - t0:.0f} s: {len(res['hits'])} hold a whole SIB1 slot, "
          f"{res['anchors']} SSB anchors -> data/{name}")
    print(f"next: matlab -batch \"addpath('matlab'); sib1_offline('data/{name}', 'FreqMHz', {a.freq}, "
          f"'SSBMHz', {a.ssb:g}, 'PCI', {a.pci})\" > logs/sib1_{name}.log 2>&1 < /dev/null")
    return 0 if res["hits"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
