"""Find NR cells: capture bursts at one or more LOs and search the whole SSB raster on the PC.

  .venv_local/bin/python -m capture.cell_scan --port /dev/ttyACM0 --freq 2131 2670
  .venv_local/bin/python -m capture.cell_scan --port /dev/ttyACM0 --freq 2g 5g          # the C5's whole range
  .venv_local/bin/python -m capture.cell_scan --port /dev/ttyACM0 --freq 2131 --ppm -5.8 --record
  .venv_local/bin/python -m capture.cell_scan --session data/<plain burst session>      # offline

Per LO: full bursts (1.64 ms at 20 MS/s), timed on the board's clock until every start time of a
20 ms-periodic SSB was held whole by --looks bursts (2: 47 bursts in the median on the board,
about 5.5 s; untimed bursts would need about 89 by simulation), then for every burst a PSS search over all GSCN positions within +-6.1 MHz of the LO
(the analog filter passes about +-8 MHz, the SSB is 3.6 MHz wide), all three N_ID2 and the CFO
hypotheses in 5 kHz steps (10 kHz at 30 kHz SCS) (nr_ssb.pss_search_raster), then SSS, PBCH and
MIB (nr_ssb.receive). Prints a table of cells and writes it as JSON.

SSB numerologies per LO (ssb_cases): 15 kHz case A below 3 GHz, plus 30 kHz case C (L_max 8)
near TDD bands (n40, n41/n38) and only that at 5 GHz (n46, n79).

The CFO is the board's crystal error times the frequency; the C5's crystal moves by several
ppm (docs/hardware.md). Without --ppm the search spans +-12 ppm (at least --cfo-span) around 0
until a decoded cell measures the crystal, then +-2 ppm around it (cfo_window).
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HALF_SPAN_HZ = 6.1e6      # +-8 MHz analog passband minus half an SSB (1.8 MHz) minus margin
SSB_PERIOD_US = 20000.0   # SSB period a UE assumes for the initial cell search (TS 38.213 4.1)
SSB_SPAN_US = 330.0       # an SSB must fit into the burst whole: 4 symbols with CP at 15 kHz (286 us) + margin
# unpaired (TDD) bands within the C5's 2 GHz range whose SSBs may use 30 kHz case C: n40, n53, n41/n38/n90
TDD_30K_MHZ = ((2300.0, 2400.0), (2483.5, 2690.0))
CFO_PPM_KNOWN = 2.0       # CFO window around the crystal's offset once a decoded cell has measured it
CFO_PPM_UNKNOWN = 12.0    # ... and around 0 before (the C5's crystal was seen from -7.3 to +2.9 ppm)


def ssb_cases(f_lo_mhz: float) -> list[tuple[float, int]]:
    """(SCS Hz, L_max) to search at an LO: 15 kHz case A (L_max 4) below 3 GHz; 30 kHz case C (L_max 8,
    TS 38.213 4.1 for unpaired spectrum above 1.88 GHz) where a TDD band lies within the LO's SSB span,
    and only that at 5 GHz (n46, n79, n96)."""
    if f_lo_mhz >= 3000:
        return [(30e3, 8)]
    half = HALF_SPAN_HZ / 1e6
    tdd = any(a - half <= f_lo_mhz <= b + half for a, b in TDD_30K_MHZ)
    return [(15e3, 4)] + ([(30e3, 8)] if tdd else [])


def cfo_window(f_mhz: float, ppm: float | None = None, span_hz: float = 30e3) -> tuple[float, float]:
    """(centre, span) of the CFO search in Hz at f_mhz: around the crystal's offset (-ppm x f) +- 2 ppm
    (at least 5 kHz) when the crystal is known, else around 0 +- 12 ppm (at least span_hz)."""
    if ppm is None:
        return 0.0, max(span_hz, CFO_PPM_UNKNOWN * f_mhz)
    return -ppm * f_mhz, max(5e3, CFO_PPM_KNOWN * f_mhz)


class Coverage:
    """The start times of an SSB, modulo its period, that the bursts so far held whole (10 us cells).
    A burst of duration d starting at t0 holds an SSB that starts in [t0, t0 + d - SSB_SPAN_US]."""

    def __init__(self, period_us: float = SSB_PERIOD_US, cell_us: float = 10.0):
        self.period, self.cell = period_us, cell_us
        self.seen = np.zeros(int(round(period_us / cell_us)), int)

    def add(self, t0_us: float, dur_us: float):
        k = int((dur_us - SSB_SPAN_US) // self.cell)
        if k > 0:
            i0 = int(np.ceil((t0_us % self.period) / self.cell))
            self.seen[(i0 + np.arange(k)) % len(self.seen)] += 1

    def looks(self) -> int:
        """How often the least seen SSB start time was held whole."""
        return int(self.seen.min())

    def next_start(self, looks: int, dur_us: float) -> float:
        """Start (us, modulo the period) that centres a burst of dur_us on the longest stretch of start
        times held fewer than `looks` times."""
        low = self.seen < looks
        if low.all() or not low.any():
            return 0.0
        shift = int(np.argmin(low))                          # a cell that is seen enough: no stretch wraps
        d = np.diff(np.concatenate([[0], np.roll(low, -shift).astype(int), [0]]))
        s, e = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
        j = int(np.argmax(e - s))
        centre = ((s[j] + e[j]) / 2 + shift) * self.cell
        return (centre - (dur_us - SSB_SPAN_US) / 2) % self.period


def sweep_eta_s(los, looks: int = 2, auto_gain: bool = True, workers: int = 2) -> float:
    """Estimated seconds for a scan of the LOs `los` with timed bursts (*estimate*, fitted to the board
    run of 24 Nov 2026: 152 LOs, 2 looks, auto gain, 2 workers, 845 s; 2 GHz alone 361 s):
    the longer of the capture (18 x looks + 11 bursts per LO, 47 for 2 looks on the board; 0.094 s
    per burst; 0.35 s per LO for tuning, 0.72 s for the gain scan) and the PC search (per burst and
    core: 0.127 s for 15 kHz and 0.273 s for 30 kHz below 3 GHz, 0.08 s for 30 kHz at 5 GHz, with the
    crystal known), plus 10 s to start."""
    bursts = 18 * looks + 11
    cpu = {(15e3, True): 0.127, (30e3, True): 0.273, (30e3, False): 0.08, (15e3, False): 0.08}
    capture = len(los) * (bursts * 0.094 + 0.35 + (0.72 if auto_gain else 0.0))
    search = sum(bursts * cpu[(scs, lo < 3000)] for lo in los for scs, _ in ssb_cases(lo)) / max(1, workers)
    return max(capture, search) + 10.0


def gscn_number(f: float) -> int:
    if f < 3000e6:
        n = int(f // 1.2e6)
        return int(3 * n + (round((f - n * 1.2e6) / 50e3) - 3) / 2)
    return int(7499 + round((f - 3000e6) / 1.44e6))


def scan_burst(x, fs: float, f_lo: float, cfo_center: float, cfo_span: float, t0_us=None, scale: float = 1.0,
               nid2s=(0, 1, 2), scs: float = 15e3, lmax: int = 4) -> list[dict]:
    """Worker: every SSB with a PCI in one burst (15 kHz case A, or 30 kHz case C). CFO hypotheses in
    steps of 5 kHz x SCS / 15 kHz (the same phase error over the shorter 30 kHz PSS symbol)."""
    from capture import nr_ssb as nr
    half = min(HALF_SPAN_HZ, fs / 2 - 1.9e6 * scs / 15e3)
    offs = [f - f_lo for f in nr.gscn_raster(f_lo - half, f_lo + half)]
    c0 = round(cfo_center / 250) * 250
    step = 5000.0 * scs / 15e3
    cfos = c0 + np.arange(-int(cfo_span // step), int(cfo_span // step) + 1) * step
    out = []
    for p in nr.pss_search_raster(x, fs, offs, cfos, nid2s=nid2s, scs=scs):
        f_ssb = f_lo + p["f_off"]
        r = nr.receive(x, fs, p["t_pss"], f_ssb, f_lo, cfo=p["cfo"], nid2=p["nid2"], t_first_us=t0_us, scale=scale,
                       scs=scs, lmax=lmax)
        if r.get("pci") is None:
            continue
        e = {"f_ssb": f_ssb, "scs": scs, "pci": r["pci"], "ok": bool(r["ok"]), "snr_db": r.get("snr_db"),
             "rsrp_dbfs": r.get("rsrp_dbfs"), "cfo": r.get("cfo"), "issb": r.get("issb"), "pss_stat": p["stat"],
             "sss_stat": r.get("sss_stat"), "t_pss_us": r.get("t_pss_us")}
        if r["ok"]:
            m = r["mib"]
            e.update(sfn=m["sfn"], mib={k: m[k] for k in ("scs_common_khz", "kssb", "dmrs_typeA_pos",
                                                          "pdcch_config_sib1", "cell_barred", "intra_freq_reselection")})
        out.append(e)
    return out


def aggregate(entries: list[dict]) -> list[dict]:
    """One row per (SSB frequency, PCI, SCS), the cells with the most MIBs first."""
    cells: dict[tuple, dict] = {}
    for e in entries:
        scs = e.get("scs", 15e3)
        k = (round(e["f_ssb"]), e["pci"], scs)
        c = cells.setdefault(k, {"f_ssb_mhz": e["f_ssb"] / 1e6, "gscn": gscn_number(e["f_ssb"]), "pci": e["pci"],
                                 "scs_khz": scs / 1e3, "ssbs": 0, "mibs": 0, "snr": [], "rsrp": [], "cfo": [],
                                 "issb": set(), "mib": None, "sfn": None})
        c["ssbs"] += 1
        c["snr"].append(e["snr_db"]); c["rsrp"].append(e["rsrp_dbfs"]); c["issb"].add(e["issb"])
        if e["ok"]:
            c["mibs"] += 1
            c["cfo"].append(e["cfo"]); c["mib"] = e["mib"]; c["sfn"] = e["sfn"]
    out = []
    for c in cells.values():
        out.append({**{k: v for k, v in c.items() if k not in ("snr", "rsrp", "cfo", "issb")},
                    "snr_db": float(np.median(c["snr"])), "rsrp_dbfs": float(np.median(c["rsrp"])),
                    "cfo_hz": float(np.median(c["cfo"])) if c["cfo"] else None, "ssb_index": sorted(c["issb"])})
    return sorted(out, key=lambda c: (-c["mibs"], -c["ssbs"]))


def capture_lo(sdr, f_lo_mhz: int, rate_code: int, n_bursts: int, pool, cfo_center: float = 0.0,
               cfo_span: float = 30e3, bits: int = 8, say=print, stop=None, writer=None, gain=None,
               cases=((15e3, 4),), looks: int | None = None) -> tuple[list, int]:
    """Capture full bursts at one LO (frequency and gain already set) and submit their scans (one per
    (SCS, L_max) case) to the pool. n_bursts bursts; or with looks, until every start time of an SSB
    within its 20 ms period was held whole by `looks` bursts (at most n_bursts). Then each burst after
    the second is delayed so that it lands on the longest stretch still missing (the spacing of the
    bursts on the board's clock varies by about +-2 ms, so it is measured as it goes).
    Returns (futures, bursts captured)."""
    from .ssb_sources import iq_scale
    cov = Coverage() if looks else None
    futs, n, gaps, last, delay, dur = [], 0, [], None, 0.0, None
    for k in range(n_bursts):
        if stop and stop():
            break
        if cov is not None and gaps and dur:
            spacing = float(np.median(gaps[-8:]))
            delay = (cov.next_start(looks, dur) - (last + spacing)) % cov.period
            time.sleep(delay / 1e6)
        else:
            delay = 0.0
        b = sdr.capture(sdr.max_samples, rate_code, bits=bits)
        if b.t0_us is None:
            cov = None                                   # firmware without CAPT: a fixed number of bursts
        elif last is not None:
            gaps.append(b.t0_us - last - delay)
        last = b.t0_us
        if not b.crc_ok:
            continue
        n += 1
        if writer:
            writer.add(b, freq_mhz=f_lo_mhz, gain=gain, bandwidth_mhz=None)
        x = b.iq()
        for scs, lmax in cases:
            futs.append(pool.submit(scan_burst, x, b.fs, f_lo_mhz * 1e6, cfo_center, cfo_span, b.t0_us,
                                    iq_scale(b.fmt), (0, 1, 2), scs, lmax))
        dur = len(x) / b.fs * 1e6
        if cov is not None:
            cov.add(b.t0_us, dur)
            if cov.looks() >= looks:
                break
        if k % 20 == 19:
            say(f"  LO {f_lo_mhz} MHz: {k + 1} bursts captured")
    return futs, n


def scan_lo(sdr, f_lo_mhz: int, rate_code: int, n_bursts: int, pool, cfo_center: float = 0.0, cfo_span: float = 30e3,
            bits: int = 8, say=print, stop=None, writer=None, gain=None, scs: float = 15e3, lmax: int = 4,
            looks: int | None = None) -> list[dict]:
    """Capture bursts at one LO (frequency and gain already set) and scan them (capture_lo)."""
    futs, _ = capture_lo(sdr, f_lo_mhz, rate_code, n_bursts, pool, cfo_center, cfo_span, bits, say, stop, writer, gain,
                         ((scs, lmax),), looks)
    return aggregate([e for f in futs for e in f.result()])


def crystal_ppm(cells: list[dict]) -> float | None:
    """The board's crystal error from the decoded cell with the most MIBs (CFO = -ppm x f), or None."""
    ok = [c for c in cells if c["mibs"] and c.get("cfo_hz") is not None]
    if not ok:
        return None
    c = max(ok, key=lambda c: c["mibs"])
    return -c["cfo_hz"] / c["f_ssb_mhz"]


def sweep(sdr, los, rate_code: int, pool, args, ppm: float | None = None, looks: int | None = 2,
          max_bursts: int = 100, say=print, stop=None, writer=None, cfo_span: float = 30e3, cfo=None, cases=None):
    """Scan LO after LO: gain (args.gain, 'auto' = pick_gain), the SSB cases of ssb_cases() (or `cases`
    for all), the CFO window of cfo_window() around the crystal measured so far (each decoded cell
    updates it; or a fixed `cfo` = (centre, span)), and bursts timed by capture_lo(looks). The PC
    searches an LO's bursts while the next LO is captured.
    Yields one dict per LO, in order: lo, gain, bursts, cases, cfo (centre, span), cells, ppm, index."""
    from .capture_bursts import _setup, pick_gain
    pending: list[dict] = []
    los = list(los)

    def finish(p):
        nonlocal ppm
        p["cells"] = aggregate([e for f in p.pop("futs") for e in f.result()])
        new = crystal_ppm(p["cells"])
        ppm = new if new is not None else ppm
        p["ppm"] = ppm
        return p

    for i, lo in enumerate(los):
        if stop and stop():
            return
        gain = args.gain
        if gain == "auto":
            gain, _ = pick_gain(sdr, lo, rate_code, 10, args)
        _setup(sdr, lo, args, int(gain))
        c0, span = cfo or cfo_window(lo, ppm, cfo_span)
        lo_cases = cases or ssb_cases(lo)
        futs, n = capture_lo(sdr, lo, rate_code, max_bursts, pool, c0, span, args.bits, say, stop, writer, int(gain),
                             lo_cases, looks)
        pending.append({"lo": lo, "gain": int(gain), "bursts": n, "cases": lo_cases, "cfo": (c0, span), "futs": futs,
                        "index": i})
        while pending and all(f.done() for f in pending[0]["futs"]):
            yield finish(pending.pop(0))
    while pending:
        if stop and stop():
            return
        yield finish(pending.pop(0))


def print_table(cells, title=""):
    print(title)
    print("  SSB (MHz)   GSCN  SCS   PCI  SSBs  MIBs  SNR dB  RSRP dBFS  CFO kHz  SSB idx  MIB")
    for c in cells:
        m = c["mib"]
        mib = (f"kSSB {m['kssb']}, pdcch {m['pdcch_config_sib1']}, " + ("barred" if m["cell_barred"] else "not barred")) if m else "-"
        cfo = f"{c['cfo_hz'] / 1e3:+7.2f}" if c["cfo_hz"] is not None else "      -"
        print(f"  {c['f_ssb_mhz']:9.3f}  {c['gscn']:5d}  {c.get('scs_khz', 15):3g}  {c['pci']:4d}  {c['ssbs']:4d}  {c['mibs']:4d}  {c['snr_db']:+6.1f}"
              f"  {c['rsrp_dbfs']:9.1f}  {cfo}  {str(c['ssb_index']):7s}  {mib}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--port")
    p.add_argument("--freq", nargs="+", help="LOs in whole MHz: 2131 2670, ranges 2175-2715/12, or 2g / 5g "
                                              "(the C5's whole 2.13-2.73 or 4.79-5.99 GHz range in 12 MHz steps)")
    p.add_argument("--session", help="scan a recorded session of full bursts instead of a board")
    p.add_argument("--looks", type=int, default=2, help="timed bursts per LO until every SSB position in the 20 ms "
                                                        "period was seen this often (0: a fixed --bursts)")
    p.add_argument("--bursts", type=int, default=100, help="bursts per LO (at most, with --looks)")
    p.add_argument("--rate", type=float, default=20, help="MS/s")
    p.add_argument("--bits", type=int, choices=(8, 10), default=8)
    p.add_argument("--gain", default="auto")
    p.add_argument("--cfo", type=float, default=None, help="fixed CFO centre (Hz) for every LO, with --cfo-span")
    p.add_argument("--ppm", type=float, default=None, help="board crystal error: CFO -ppm x f +- 2 ppm per LO")
    p.add_argument("--cfo-span", type=float, default=30e3, help="CFO span (Hz) with --cfo, and the least one before "
                                                                 "the crystal is known")
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--scs", type=float, default=None, choices=(15, 30),
                   help="SSB subcarrier spacing in kHz for every LO (default: per LO, cell_scan.ssb_cases)")
    p.add_argument("--lmax", type=int, default=None, choices=(4, 8), help="SSBs per half frame (default 4 at 15 kHz, 8 at 30)")
    p.add_argument("--record", action="store_true", help="store the bursts as a session in data/")
    p.add_argument("--out", default=None, help="JSON with the cells (default results/cells_<time>.json)")
    a = p.parse_args(argv)
    if not a.session and not (a.port and a.freq):
        p.error("give --port and --freq, or --session")
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    from .ssb_sources import _worker_init, iq_scale, parse_los
    forced = None
    if a.scs:
        forced = [(a.scs * 1e3, a.lmax or (4 if a.scs == 15 else 8))]
    pool = ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context("spawn"), initializer=_worker_init)
    result = {"created": time.strftime("%Y-%m-%dT%H:%M:%S"), "los": []}
    t0 = time.monotonic()
    try:
        if a.session:
            from .session import load_session
            meta, recs, iq = load_session(a.session)
            by_lo: dict[float, list] = {}
            for r, x in zip(recs, iq):
                lo = float(r["freq_mhz"])
                c, span = (a.cfo, a.cfo_span) if a.cfo is not None else cfo_window(lo, a.ppm, a.cfo_span)
                for scs, lmax in forced or ssb_cases(lo):
                    by_lo.setdefault(lo, []).append(pool.submit(scan_burst, x, r["fs"], lo * 1e6, c, span, r.get("t0_us"),
                                                                iq_scale(r["fmt"]), (0, 1, 2), scs, lmax))
            for lo, futs in sorted(by_lo.items()):
                cells = aggregate([e for f in futs for e in f.result()])
                print_table(cells, f"LO {lo:g} MHz, {len(futs)} burst searches:")
                result["los"].append({"lo_mhz": lo, "bursts": len(futs), "cells": cells})
        else:
            from .esp_sdr import EspSdr, rate_code_for
            from .session import SessionWriter
            los = parse_los(" ".join(a.freq))
            sdr = EspSdr(a.port, timeout=3.0)
            device = sdr.connect()
            rate = rate_code_for(a.rate * 1e6)
            a.bandwidth, a.dccal = None, False
            writer = None
            if a.record:
                writer = SessionWriter(ROOT / "data" / f"{time.strftime('%Y%m%d_%H%M%S')}_cellscan", device,
                                       {"tool": "cell_scan", "freq_mhz": los, "rate_msps": a.rate, "bits": a.bits,
                                        "gain": a.gain, "looks": a.looks})
            print(f"{len(los)} LOs, about {sweep_eta_s(los, a.looks or 4, a.gain == 'auto', a.workers) / 60:.0f} min "
                  f"(estimate)")
            try:
                for r in sweep(sdr, los, rate, pool, a, ppm=a.ppm, looks=a.looks or None, max_bursts=a.bursts,
                               writer=writer, cfo_span=a.cfo_span, cfo=None if a.cfo is None else (a.cfo, a.cfo_span),
                               cases=forced):
                    print_table(r["cells"], f"LO {r['lo']} MHz (gain {r['gain']}), {r['bursts']} bursts, "
                                            f"CFO {r['cfo'][0] / 1e3:+.1f} +- {r['cfo'][1] / 1e3:.1f} kHz, "
                                            f"SCS {'/'.join(f'{s / 1e3:g}' for s, _ in r['cases'])} kHz:")
                    result["los"].append({"lo_mhz": r["lo"], "gain": r["gain"], "bursts": r["bursts"], "cells": r["cells"],
                                          "cfo_hz": list(r["cfo"]), "ppm": r["ppm"],
                                          "scs_khz": [s / 1e3 for s, _ in r["cases"]]})
            finally:
                if writer:
                    writer.close()
                sdr.close()
    finally:
        pool.shutdown()
    out = Path(a.out) if a.out else ROOT / "results" / f"cells_{time.strftime('%Y%m%d_%H%M%S')}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1))
    print(f"{sum(len(x['cells']) for x in result['los'])} cells in {time.monotonic() - t0:.0f} s -> {out}")

if __name__ == "__main__":
    main()
