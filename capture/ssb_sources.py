"""Board, replay and decoding for the SSB tools (no GUI): used by the live view (capture/ssb_live.py),
the headless monitor (capture/ssb_monitor.py) and the cell scanner (capture/cell_scan.py).

- Channel helpers: the C5's tuning ranges, NR bands, LO lists (parse_los), an LO for an SSB.
- Decoder worker: decode_job runs nr_ssb.receive on one cut in a separate process.
- Sources: threads that put events (frames, logs, scan results) on a queue: LiveSource (the board
  with SSBDET, reconnects by itself), ScanSource (a scan over LOs from the window), ReplaySource
  (a recorded session).
- Command-line options shared by the live view and the monitor (add_source_args, check_source_args).

Receive only: the board never transmits.
"""
from __future__ import annotations

import argparse
import os
import queue
import threading
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MAX_BACKLOG = 16           # decode jobs in flight before live cuts are skipped (covers one redraw)
EXTRA_SSS_THRESHOLD = 18.0  # other N_ID2 tried on every cut: kept if the MIB passes the CRC, or else
                           # only above this SSS statistic (336 e^-18 = 5e-6 false PCIs per try; at the
                           # normal 13 about one false PCI per 800 cuts: PCI 647 on 23 Nov)
C5_TUNE_MHZ = ((2130, 2730), (4790, 5990))   # LOs the C5 tunes exactly (docs/hardware.md)
SSB_HALF_SPAN_MHZ = 6.1                      # SSB centre at most this far from the LO (cell_scan.HALF_SPAN_HZ)
SCAN_BURSTS = 100                            # full bursts per LO at most for a scan from the window
SCAN_LOOKS = 2                               # ... timed until every SSB position of a 20 ms period was seen twice (cell_scan.capture_lo)
SCAN_LO_STEP = 12                            # MHz between scan LOs: SSB centres +-6.1 MHz around each, without gaps
SWEEPS = {"2g": C5_TUNE_MHZ[0], "5g": C5_TUNE_MHZ[1]}   # the C5's whole ranges, for parse_los
# NR downlink ranges inside the C5's tuning ranges (TS 38.101-1); n7 and n41 overlap
NR_BANDS = (("n1", 2110, 2170), ("n40", 2300, 2400), ("n41", 2496, 2690), ("n38", 2570, 2620), ("n7", 2620, 2690),
            ("n79", 4400, 5000), ("n46", 5150, 5925), ("n47", 5855, 5925))


def nr_band(f_mhz: float) -> str:
    """NR band name(s) of a downlink frequency, e.g. 'n1' or 'n41/n7'."""
    return "/".join(b for b, lo, hi in NR_BANDS if lo <= f_mhz <= hi) or "-"


def lo_tunable(lo_mhz: int) -> bool:
    return any(a <= lo_mhz <= b for a, b in C5_TUNE_MHZ)


def parse_los(text: str) -> list[int]:
    """LOs in whole MHz from '2131 2670', '2131, 2670', ranges '2175-2715/15' (step 12 if omitted), or
    '2g' / '5g' for the C5's whole 2130-2730 / 4790-5990 MHz in 12 MHz steps."""
    out: list[int] = []
    for tok in text.replace(",", " ").lower().split():
        if tok in SWEEPS:
            a, b = SWEEPS[tok]
            out += list(range(a, b + 1, SCAN_LO_STEP))
        elif "-" in tok:
            rng, _, step = tok.partition("/")
            a, b = (int(float(v)) for v in rng.split("-", 1))
            out += list(range(a, b + 1, int(step or SCAN_LO_STEP)))
        else:
            out.append(int(float(tok)))
    return out


def lo_for_ssb(ssb_mhz: float, lo_mhz: int | None = None) -> int:
    """An LO (whole MHz) for an SSB centre: the given one if the SSB lies within +-6.1 MHz and away
    from the LO's DC (2 MHz), else about 5 MHz above the SSB (below near the top of a range)."""
    if lo_mhz is not None and 2.0 <= abs(ssb_mhz - lo_mhz) <= SSB_HALF_SPAN_MHZ and lo_tunable(lo_mhz):
        return int(lo_mhz)
    for lo in (round(ssb_mhz) + 5, round(ssb_mhz) - 5, round(ssb_mhz) + 4, round(ssb_mhz) - 4):
        if lo_tunable(lo):
            return int(lo)
    raise ValueError(f"no tunable LO within {SSB_HALF_SPAN_MHZ} MHz of {ssb_mhz} MHz (C5: {C5_TUNE_MHZ})")


def iq_scale(fmt: str, R: int = 1, gsh: int = 0) -> float:
    """Factor from payload units to 10-bit ADC counts (full scale 512), for RSRP in dBFS:
    CAP16 carries the top 8 of 10 bits (x4); a decimated cut (ssb3) is the DDC output, 3 bits
    above the input (/8), shifted right by gsh (x 2^gsh)."""
    s = 4.0 if fmt == "CAP16" else 1.0
    return s * (2.0 ** gsh) / 8.0 if R > 1 else s


# ============================================================================ worker process
def _worker_init():
    try:
        os.nice(5)
    except OSError:
        pass


def decode_job(job: dict) -> list[dict]:
    """Runs in the worker process: PSS search if needed, then nr_ssb.receive per SSB."""
    from capture import nr_ssb as nr
    t0 = time.perf_counter()
    x = job["iq"]
    fs = job["fs"]
    if job.get("t_pss") is None:
        cands = [(p["t_pss"], p["nid2"], p["cfo"], p["stat"])
                 for p in nr.pss_search(x, fs, job["f_ssb"] - job["f_lo"], nid2s=job.get("nid2s", (0, 1, 2)))]
    else:
        # the chip's cell first; then the other two N_ID2 at the same time: sectors of the same
        # site (PCI 166 next to 165) send their SSBs together. The SSS threshold keeps out noise.
        cands = [(job["t_pss"], k, job["cfo"] or 0.0, job.get("stat") if k == job["nid2"] else None)
                 for k in [job["nid2"]] + [k for k in (0, 1, 2) if k != job["nid2"]]]
    # ssb3 decimated cut: at fs / R already, SSB centre at 0 Hz, cfo_ddc removed by the chip
    cfo_ddc = job.get("cfo_ddc") or 0
    f_lo = job["f_ssb"] if (job.get("R") or 1) > 1 else job["f_lo"]
    out = []
    for t_pss, nid2, cfo, stat in cands:
        r = nr.receive(x, fs, t_pss, job["f_ssb"], f_lo, cfo=cfo - cfo_ddc, nid2=nid2, t_first_us=job.get("t_first_us"),
                       scale=job.get("scale") or 1.0)
        if r.get("cfo") is not None:
            r["cfo"] += cfo_ddc
        if stat is None and job.get("t_pss") is not None and (r.get("pci") is None or (
                not r["ok"] and (r.get("sss_stat") or 0) < EXTRA_SSS_THRESHOLD)):
            continue                                   # other N_ID2: nothing there
        keep = {k: r.get(k) for k in ("ok", "pci", "nid2", "issb", "hrf_dmrs", "cfo", "snr_db", "sss_stat", "mib",
                                      "t_pss_us", "t_frame_us", "grid", "eq", "channel", "nvar", "error", "t_pss",
                                      "rsrp_dbfs")}
        keep["pci_chip"] = job.get("pci_chip") if stat is not None else None
        keep.update(seq=job.get("seq"), mode=job.get("mode"), stat=stat, cfo_in=cfo, f_ssb=job["f_ssb"],
                    epoch=job.get("epoch", 0), run=job.get("run"))
        if r["ok"]:
            # after decoding, the sent PBCH symbols are known: re-encode the MIB (as nrBCH + nrPBCH)
            m = r["mib"]
            sent = nr.bch_encode(np.array([int(b) for b in m["bits"]], np.uint8), m["sfn"], m["hrf"], m["kssb"],
                                 r["pci"]) ^ nr.pbch_scrambling(r["pci"], r["issb"])
            keep["sent_idx"] = 2 * sent[0::2] + sent[1::2]          # 0: +,+  1: +,-  2: -,+  3: -,-
            eq = r["eq"]
            hard = np.empty(864, np.uint8)
            hard[0::2], hard[1::2] = eq.real < 0, eq.imag < 0
            keep["raw_agree"] = float(np.mean(hard == sent))
        # power envelope of the cut around the SSB (0.5 us averages), time from the PSS useful start
        tp = r.get("t_pss", t_pss)
        a = np.abs(x - np.mean(x)) ** 2
        w = max(1, int(round(0.5e-6 * fs)))
        i0 = max(0, int(tp - 30e-6 * fs))
        i1 = min(len(a), int(tp + 320e-6 * fs))
        seg = a[i0:i1][: (i1 - i0) // w * w].reshape(-1, w).mean(axis=1)
        keep["env_t"] = (i0 + w * np.arange(len(seg)) + w / 2 - tp) / fs * 1e6
        keep["env_db"] = 10 * np.log10(np.maximum(seg, 1e-3))
        out.append(keep)
    for k in out:
        k["decode_ms"] = 1e3 * (time.perf_counter() - t0) / max(1, len(out))
    return out


def scan_burst(x: np.ndarray, fs: float, f_ssb: float, f_lo: float, nid2s) -> list[tuple]:
    """Worker: (pci, pss statistic, cfo, mib ok) of every SSB found in one full burst."""
    from capture import nr_ssb as nr
    found = []
    for p in nr.pss_search(x, fs, f_ssb - f_lo, nid2s=nid2s):
        r = nr.receive(x, fs, p["t_pss"], f_ssb, f_lo, cfo=p["cfo"], nid2=p["nid2"])
        if r.get("pci") is not None:
            found.append((r["pci"], p["stat"], r["cfo"], r["ok"]))
    return found


# ============================================================================ sources
class Source(threading.Thread):
    def __init__(self, args, out: queue.Queue, pool: ProcessPoolExecutor):
        super().__init__(daemon=True)
        self.args, self.out, self.pool = args, out, pool
        self.stop_flag = threading.Event()

    def say(self, text: str):
        self.out.put({"kind": "log", "text": text})

    def run(self):
        try:
            self.work()
        except Exception as e:  # shown in the window, not lost on stderr
            self.out.put({"kind": "error", "text": f"{type(e).__name__}: {e}"})
        finally:
            self.out.put({"kind": "stopped"})


class LiveSource(Source):
    """The board with SSBDET. If the board drops out (USB re-enumeration after a reset of its
    own) or stops answering, it is reopened, reset through RTS if needed, and SSBDET restarts
    with the last measured CFO; the session (if recorded) continues."""

    def work(self):
        from .capture_bursts import _setup, pick_gain
        from .esp_sdr import EspSdrError, rate_code_for, ssb_stream
        from .session import SsbRecorder
        a = self.args
        rate = rate_code_for(a.rate * 1e6)
        offset = None if a.ssb == "auto" else int(round((a.ssb - a.freq) * 1e6))
        gain, nid2, cfo = a.gain, a.nid2, a.cfo
        rec, configured, idle_runs = None, False, 0
        try:
            while not self.stop_flag.is_set():
                sdr, device = self.open_board()
                if sdr is None:
                    return
                try:
                    if device.get("boot"):
                        self.say(f"board start reason: {device['boot']}")      # e.g. brownout after a reset
                    v3 = "SSBDET3" in sdr.caps
                    if not configured:
                        self.say(f"{sdr.info_line}  {device.get('fwopt') or ''}")
                        if "SSBDET" not in sdr.caps:
                            raise RuntimeError("this firmware has no SSBDET: flash the ssb-sniffer build (firmware/fw.sh flash)")
                        if gain == "auto":
                            self.say("auto gain: scanning gain indices")
                            gain, _ = pick_gain(sdr, a.freq, rate, 10, a)
                        gain = int(gain)
                    _setup(sdr, a.freq, a, gain)
                    if offset is None:
                        # --ssb auto: search the whole SSB raster around the LO on the PC first
                        from .cell_scan import scan_lo
                        self.say(f"cell scan: all GSCN positions within +-6.1 MHz of {a.freq} MHz, "
                                 f"{a.scan_bursts} bursts, CFO {cfo or 0:+.0f} +- {a.cfo_span} Hz")
                        cells = scan_lo(sdr, a.freq, rate, a.scan_bursts, self.pool, cfo or 0, a.cfo_span,
                                        say=self.say, stop=self.stop_flag.is_set, looks=SCAN_LOOKS)
                        if self.stop_flag.is_set():
                            return
                        if not cells:
                            raise RuntimeError("cell scan found no SSB near this LO: try another --freq")
                        for c in cells[:5]:
                            self.say(f"  {c['f_ssb_mhz']:.3f} MHz PCI {c['pci']}: {c['ssbs']} SSBs, {c['mibs']} MIBs, "
                                     f"SNR {c['snr_db']:+.1f} dB")
                        best = cells[0]
                        a.ssb = best["f_ssb_mhz"]
                        offset = int(round((a.ssb - a.freq) * 1e6))
                        nid2, cfo = best["pci"] % 3, best["cfo_hz"] if best["cfo_hz"] is not None else (cfo or 0)
                        self.say(f"tracking PCI {best['pci']} at {a.ssb:.3f} MHz")
                    if not configured:
                        self.say(f"LO {a.freq} MHz, gain index {gain}")
                        if (nid2 is None or cfo is None) and not v3:
                            nid2, cfo = self.scan(sdr, rate, nid2)
                            if self.stop_flag.is_set():
                                return
                        cfg = {"kind": "config", "source": f"live {a.port}", "f_lo": a.freq * 1e6, "f_ssb": a.ssb * 1e6,
                               "rate": a.rate, "bits": a.bits, "gain": gain, "nid2": nid2, "cfo": cfo, "record": None}
                        if a.record:
                            stamp = time.strftime("%Y%m%d_%H%M%S")
                            settings = {"mode": "ssbdet", "tool": "ssb_live", "freq_mhz": [a.freq], "ssb_mhz": a.ssb,
                                        "offset_hz": offset, "nid2": nid2, "cfo_hz": cfo, "rate_msps": a.rate,
                                        "rate_code": rate, "bits": a.bits, "gain": gain, "bandwidth_mhz": a.bandwidth,
                                        "thr": a.thr, "pre_us": a.pre, "post_us": a.post, "period_ms": a.period,
                                        "every": a.every, "track": True, "label": a.label, "note": a.note}
                            rec = SsbRecorder(ROOT / "data" / f"{stamp}_{a.label or 'ssblive'}", device, settings,
                                              a.freq, gain, a.bandwidth, rate)
                            cfg["record"] = str(rec.folder.relative_to(ROOT))
                        self.out.put(cfg)
                        configured = True
                    elif idle_runs >= 2 and not v3:
                        self.say("no SSB in two runs: scanning for the cell again (the CFO may have drifted)")
                        nid2, cfo = self.scan(sdr, rate, None)
                        idle_runs = 0
                        if self.stop_flag.is_set():
                            return
                    idle = False
                    # ssb3: the chip finds the cell itself (all N_ID2 when unknown, CFO +- span: the crystal
                    # moves by several ppm across resets), and sends decimated cuts (4x less data). A cell
                    # picked from a scan comes with its CFO: a narrower span (track_span) finds it sooner.
                    span = (getattr(a, "track_span", None) or a.cfo_span) if v3 else 0
                    run_nid2 = (3 if nid2 is None else nid2) if v3 else nid2
                    if v3 and nid2 is None:
                        self.say(f"cell search on the chip: all N_ID2, CFO {cfo or 0:+.0f} +- {span} Hz")
                    for kind, item in ssb_stream(sdr, rate, offset, run_nid2, int(round(cfo or 0)), int(round(a.thr * 10)),
                                                 a.pre, a.post, int(a.period * 1000), a.bits, a.every, 0, True, None,
                                                 stop=self.stop_flag.is_set, idle_timeout=None if v3 else 15.0,
                                                 dec=v3 and not a.raw_cuts, cfo_span=span):
                        if kind == "idle":
                            idle = True
                            idle_runs += 1
                            self.say(f"no SSB for 15 s (board alive): restarting SSBDET ({item})")
                            continue
                        if kind == "text":
                            self.say(f"board: {item}")
                            continue
                        if kind != "frame":
                            self.say(item)
                            continue
                        f = item
                        idle_runs = 0
                        cfo = f.cfo_hz                 # restart from here after a reset: the crystal drifts
                        nid2 = f.nid2
                        if rec:
                            rec.add(f)
                        self.out.put({"kind": "frame", "seq": f.seq, "mode": f.mode, "stat": f.stat, "nid2": f.nid2,
                                      "cfo": f.cfo_hz, "fs": f.fs_cut, "crc_ok": f.crc_ok, "host_time": f.host_time,
                                      "t_pss_chip_us": f.t0_us + f.t0_sub_ns / 1000 + (f.flip + f.pss_idx) / f.fs * 1e6,
                                      "iq": f.iq() if f.n and f.crc_ok else None, "t_pss": f.pss_in_cut,
                                      "t_first_us": f.t_first_us, "f_ssb": a.ssb * 1e6, "f_lo": a.freq * 1e6,
                                      "R": f.R, "cfo_ddc": f.cfo_ddc_hz, "pci_chip": f.pci,
                                      "scale": iq_scale(f.fmt, f.R, f.gsh)})
                    if not idle:
                        return                                   # stopped normally
                except (OSError, TimeoutError, EspSdrError) as e:  # serial.SerialException is an OSError
                    if self.stop_flag.is_set():
                        return
                    self.out.put({"kind": "reconnect"})
                    self.say(f"board lost ({type(e).__name__}: {str(e)[:80]}); reconnecting")
                finally:
                    try:
                        sdr.close()
                    except Exception:
                        pass
        finally:
            if rec:
                rec.close()
                self.say(f"recorded {rec.count} cuts -> {rec.folder.relative_to(ROOT)}")

    def open_board(self):
        """Open and SYNC the board; if it does not answer, reset it through RTS and try again.
        Waits for the port to come back after a USB re-enumeration. Returns (sdr, device) or (None, None)."""
        from .esp_sdr import EspSdr, EspSdrError
        a = self.args
        for attempt in range(1, 7):
            if self.stop_flag.is_set():
                return None, None
            t_end = time.monotonic() + 15
            while not Path(a.port).exists() and time.monotonic() < t_end and not self.stop_flag.is_set():
                time.sleep(0.2)
            if not Path(a.port).exists():
                # after a USB re-enumeration the board can come back as another ttyACM: use its fixed id
                ids = sorted(Path("/dev/serial/by-id").glob("usb-Espressif*"))
                if ids:
                    self.say(f"{a.port} is gone; using {ids[0].name} -> {ids[0].resolve().name}")
                    a.port = str(ids[0])
            self.say(f"connecting to {a.port}" + (f" (attempt {attempt})" if attempt > 1 else ""))
            try:
                sdr = EspSdr(a.port, timeout=3.0)
            except OSError as e:
                self.say(f"cannot open {a.port}: {e}")
                time.sleep(1.5)
                continue
            for reset in (False, True):
                try:
                    if reset:
                        self.say("board does not answer: resetting it (RTS, nothing is written)")
                        sdr.hard_reset()
                    return sdr, sdr.connect()
                except (OSError, TimeoutError, EspSdrError):
                    continue
            try:
                sdr.ser.close()
            except Exception:
                pass
            time.sleep(1.5)
        raise RuntimeError(f"board on {a.port} not reachable after 6 attempts (unplug and replug it?)")

    def scan(self, sdr, rate, nid2):
        """Full bursts at the SSB frequency, searched on the PC for all N_ID2 (or the given one)
        and +-30 kHz: returns (N_ID2, CFO) of the cell with the most decoded MIBs."""
        a = self.args
        nid2s = (0, 1, 2) if nid2 is None else (nid2,)
        self.say(f"scan: up to {a.scan_bursts} full bursts, N_ID2 {list(nid2s)}, CFO -30..+30 kHz")
        futs, found = [], {}
        for k in range(a.scan_bursts):
            if self.stop_flag.is_set():
                return None, None
            b = sdr.capture(sdr.max_samples, rate, bits=8)
            if b.crc_ok:
                futs.append(self.pool.submit(scan_burst, b.iq(), b.fs, a.ssb * 1e6, a.freq * 1e6, nid2s))
            done = [f for f in futs if f.done()]
            for f in done:
                futs.remove(f)
                for pci, stat, cfo, ok in f.result():
                    found.setdefault(pci, []).append((stat, cfo, ok))
            mibs = {p: sum(v[2] for v in found[p]) for p in found}
            if mibs and max(mibs.values()) >= 3:
                break
            if k % 10 == 9:
                self.say(f"scan: {k + 1} bursts, " + (", ".join(f"PCI {p}: {len(v)} SSBs, {mibs[p]} MIBs"
                                                                  for p, v in found.items()) or "nothing yet"))
        for f in futs:
            for pci, stat, cfo, ok in f.result():
                found.setdefault(pci, []).append((stat, cfo, ok))
        if not found:
            raise RuntimeError("scan found no SSB at this frequency; check --ssb, or give --nid2 and --cfo")
        pci = max(found, key=lambda p: (sum(v[2] for v in found[p]), len(found[p])))
        cfo = float(np.median([v[1] for v in found[pci]]))
        self.say(f"scan: strongest cell PCI {pci} ({len(found[pci])} SSBs, {sum(v[2] for v in found[pci])} MIBs), "
                 f"CFO {cfo:+.0f} Hz -> SSBDET with N_ID2 {pci % 3}")
        return pci % 3, cfo


class ScanSource(LiveSource):
    """Cell scan started from the window (cell_scan.sweep): at each LO, full bursts timed until every
    SSB position of a 20 ms period was seen `looks` times, searched on the PC over the whole SSB raster
    (+-6.1 MHz), all N_ID2, the SSB numerologies of the band, and a CFO window around the crystal
    (narrowed once a decoded cell has measured it). Puts one {"kind": "scan"} event per LO with its
    cells and the crystal. Shares the board handling of LiveSource."""

    def __init__(self, args, out, pool, los, ppm=None, looks: int | None = SCAN_LOOKS, n_bursts: int = SCAN_BURSTS):
        super().__init__(args, out, pool)
        self.los, self.ppm, self.looks, self.n_bursts = list(los), ppm, looks, n_bursts

    def work(self):
        from .cell_scan import sweep
        from .esp_sdr import rate_code_for
        a = self.args
        rate = rate_code_for(a.rate * 1e6)
        sdr, _ = self.open_board()
        if sdr is None:
            return
        t0, n = time.monotonic(), len(self.los)
        try:
            for r in sweep(sdr, self.los, rate, self.pool, a, ppm=self.ppm, looks=self.looks, max_bursts=self.n_bursts,
                           say=lambda _t: None, stop=self.stop_flag.is_set, cfo_span=a.cfo_span):
                done = r["index"] + 1
                left = (time.monotonic() - t0) / done * (n - done)
                self.out.put({"kind": "scan", **r, "n": n, "left_s": left})
        finally:
            try:
                sdr.close()
            except Exception:
                pass


class ReplaySource(Source):
    def __init__(self, args, out, pool, busy):
        super().__init__(args, out, pool)
        self.busy = busy

    def work(self):
        from .session import load_session
        a = self.args
        folder = Path(a.replay)
        meta, recs, iq = load_session(folder)
        st = meta.get("settings", {})
        f_lo = float((st.get("freq_mhz") or [0])[0]) * 1e6
        ssb = a.ssb or st.get("ssb_mhz")
        if not ssb:
            raise RuntimeError("the session does not name the SSB frequency: give --ssb")
        f_ssb = float(ssb) * 1e6
        ssbdet = bool(recs) and "ssb_pss_in_cut" in recs[0]
        self.out.put({"kind": "config", "source": f"replay {folder.name}", "f_lo": f_lo, "f_ssb": f_ssb,
                      "rate": recs[0]["fs"] / 1e6 if recs else None, "bits": 8 if recs and recs[0]["fmt"] == "CAP16" else 10,
                      "gain": st.get("gain"), "nid2": st.get("nid2"), "cfo": st.get("cfo_hz"), "record": None})
        self.say(f"replay of {len(recs)} {'SSBDET cuts' if ssbdet else 'bursts (PSS search on the PC)'}, "
                 f"speed {a.speed:g}x" + (" (as fast as decoded)" if a.speed == 0 else ""))
        wall0, t_first = time.monotonic(), None
        for rec, x in zip(recs, iq):
            if self.stop_flag.is_set():
                break
            t = rec.get("t0_us")
            if a.speed > 0 and t is not None:
                t_first = t if t_first is None else t_first
                target = wall0 + (t - t_first) / 1e6 / a.speed
                while time.monotonic() < target and not self.stop_flag.is_set():
                    time.sleep(min(0.05, target - time.monotonic()))
            while self.busy() and not self.stop_flag.is_set():   # replay never skips cuts
                time.sleep(0.02)
            R = rec.get("ssb_R", 1) or 1
            ev = {"kind": "frame", "seq": rec.get("seq", rec["idx"]), "mode": rec.get("ssb_mode", "R"),
                  "stat": rec.get("ssb_stat"), "nid2": rec.get("ssb_nid2"), "cfo": rec.get("ssb_cfo_hz"),
                  "fs": rec["fs"], "crc_ok": True, "host_time": time.time(), "iq": x,
                  "t_pss": rec.get("ssb_pss_in_cut") if ssbdet else None, "t_first_us": t,
                  "t_pss_chip_us": (t + rec["ssb_pss_in_cut"] / rec["fs"] * 1e6) if ssbdet and t is not None else None,
                  "f_ssb": f_ssb, "f_lo": f_lo, "R": R, "cfo_ddc": rec.get("ssb_cfo_ddc_hz", 0),
                  "pci_chip": rec.get("ssb_pci_chip"), "scale": iq_scale(rec["fmt"], R, rec.get("ssb_gsh", 0))}
            self.out.put(ev)
        while self.busy() and not self.stop_flag.is_set():
            time.sleep(0.05)
        self.say("replay finished")


# ============================================================================ command line
def add_source_args(p: argparse.ArgumentParser, label: str = "ssblive"):
    """Board / replay options shared by the live view and the headless monitor."""
    p.add_argument("--port", help="board, e.g. /dev/ttyACM0 (live mode)")
    p.add_argument("--replay", help="session folder to replay instead of a board")
    p.add_argument("--speed", type=float, default=1.0, help="replay speed (1 = real time, 0 = as fast as decoded)")
    p.add_argument("--freq", type=int, help="LO in whole MHz (live)")
    p.add_argument("--ssb", help="SSB centre (SS_REF) in MHz, e.g. 2125.45, or 'auto' (scan the raster first)")
    p.add_argument("--nid2", type=int, default=None, help="PSS root of the cell; omit to search")
    p.add_argument("--cfo", type=float, default=None, help="CFO in Hz to start from (crystal); omit to search")
    p.add_argument("--gain", default="auto", help="manual gain index or 'auto'")
    p.add_argument("--rate", type=float, default=20, help="MS/s")
    p.add_argument("--bits", type=int, choices=(8, 10), default=8)
    p.add_argument("--thr", type=float, default=14.0, help="on-chip PSS threshold (noise: mean 1)")
    p.add_argument("--pre", type=int, default=20)
    p.add_argument("--post", type=int, default=20)
    p.add_argument("--period", type=float, default=20.0, help="SSB period of the cell in ms")
    p.add_argument("--every", type=int, default=1, help="I/Q cut of every k-th SSB")
    p.add_argument("--scan-bursts", type=int, default=150, help="full bursts at most for the PC scan (firmware before ssb3)")
    p.add_argument("--cfo-span", type=int, default=30000, help="ssb3: CFO search span of the chip's cell search (Hz)")
    p.add_argument("--raw-cuts", action="store_true", help="ssb3: raw cuts instead of decimated 5 MS/s cuts")
    p.add_argument("--record", action="store_true", help="store the cuts as a session in data/")
    p.add_argument("--label", default=label)
    p.add_argument("--note", default="")
    p.add_argument("--workers", type=int, default=2, help="decoder processes")


def check_source_args(p: argparse.ArgumentParser, a):
    if not a.replay and not (a.port and a.freq):
        p.error("live mode needs --port and --freq (or use --replay)")
    if a.ssb is None and not a.replay:
        a.ssb = "auto"
    if a.ssb not in (None, "auto"):
        a.ssb = float(a.ssb)
    if a.replay and a.ssb == "auto":
        p.error("--ssb auto needs a board")
    a.bandwidth, a.dccal = None, False                  # used by capture_bursts._setup
