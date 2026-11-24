"""Live view of NR SS/PBCH blocks: ESP32-C5 with on-chip SSB detection (SSBDET), or a replay.

Run from the project root:

  .venv_local/bin/python -m capture.ssb_live --port /dev/ttyACM0 --freq 2131 --ssb 2125.45
  .venv_local/bin/python -m capture.ssb_live --port /dev/ttyACM0 --freq 2131 --ssb 2125.45 --nid2 0 --cfo 15000 --record
  .venv_local/bin/python -m capture.ssb_live --replay data/20261123_140216_ssbdet_v2
  .venv_local/bin/python -m capture.ssb_live --replay data/20261123_112009_n1_2131 --ssb 2125.45

Live: the board (firmware branch ssb-sniffer) finds the SSBs of one cell and sends the raw I/Q
around each; this tool decodes every cut on the PC (capture/nr_ssb.py: PCI, SSB index, SFN,
MIB, CFO, SNR, channel) and shows it. Without --nid2/--cfo it first scans with full bursts
at the SSB frequency and takes the cell with the most decoded MIBs. --record stores the cuts
as a session (same format as capture_bursts ssb-stream).
Replay: an SSBDET session (cuts with the chip's PSS positions) or a session of plain bursts
(searched for PSS on the PC; give --ssb if the session does not name the SSB frequency).

Decoding runs in a separate process. When it falls behind in live mode, cuts are counted
as skipped instead of queued, so the view stays current.
Receive only: the board never transmits.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import queue
import threading
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MAX_BACKLOG = 16           # decode jobs in flight before live cuts are skipped (covers one redraw)
EXTRA_SSS_THRESHOLD = 18.0  # other N_ID2 tried on every cut: kept if the MIB passes the CRC, or else
                           # only above this SSS statistic (336 e^-18 = 5e-6 false PCIs per try; at the
                           # normal 13 about one false PCI per 800 cuts: PCI 647 on 23 Nov)
REDRAW_MS = 700            # a redraw takes 130-300 ms of the GUI process; results queue up meanwhile
PCI_COLORS = ("#2a6fdb", "#e07b00", "#1a9e5a", "#b23a9a", "#7a6a00", "#5a5a5a")
SYM_COLORS = ("#2a6fdb", "#e07b00", "#1a9e5a", "#b23a9a")    # sent QPSK symbol (+,+) (+,-) (-,+) (-,-)


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
                    epoch=job.get("epoch", 0))
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
                                        say=self.say, stop=self.stop_flag.is_set)
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
                    # moves by several ppm across resets), and sends decimated cuts (4x less data)
                    span = a.cfo_span if v3 else 0
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


# ============================================================================ the window
class App:
    def __init__(self, root, args):
        import tkinter as tk
        from tkinter import ttk

        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure

        self.tk, self.ttk = tk, ttk
        self.root, self.args = root, args
        self.events: queue.Queue = queue.Queue()
        # one BLAS thread per decoder process: the default (a thread per core, spinning while
        # idle) made two workers take 200 ms per cut instead of 7 ms (23 Nov)
        for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ[var] = "1"
        self.pool = ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn"),
                                        initializer=_worker_init)
        self.results: queue.Queue = queue.Queue()    # filled by the pool's callbacks
        self.inflight = 0
        self.lock = threading.Lock()
        self.source: Source | None = None
        self.draw_ms = float("nan")
        self.paused = tk.BooleanVar(value=False)
        self.hidden: set[int] = set(args.hide or [])  # PCIs not shown (kept over Start/Stop)
        self._skip_release = False
        self.dirty = False
        self.reset_data()

        root.title("SSB live view (ESP32-C5)")
        root.geometry(args.geometry)
        root.columnconfigure(0, weight=1)
        root.rowconfigure(1, weight=1)
        bar = ttk.Frame(root, padding=(8, 6))
        bar.grid(row=0, column=0, columnspan=2, sticky="ew")
        self.btn_start = ttk.Button(bar, text="Start", command=self.start)
        self.btn_stop = ttk.Button(bar, text="Stop", command=self.stop, state="disabled")
        self.btn_start.pack(side="left")
        self.btn_stop.pack(side="left", padx=(4, 12))
        ttk.Checkbutton(bar, text="Pause plots", variable=self.paused).pack(side="left")
        ttk.Button(bar, text="Save PNG", command=self.save_png).pack(side="left", padx=(12, 0))
        self.src_text = tk.StringVar(value=self.describe_args())
        ttk.Label(bar, textvariable=self.src_text).pack(side="left", padx=16)

        self.fig = Figure(figsize=(11, 8), dpi=args.dpi, layout="constrained")
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        self.build_axes()

        side = ttk.Frame(root, padding=(6, 4), width=450)
        side.grid(row=1, column=1, sticky="ns")
        side.grid_propagate(False)
        side.rowconfigure(3, weight=1)
        mono = ("DejaVu Sans Mono", 9)
        st = ttk.LabelFrame(side, text="Stream", padding=4)
        st.grid(row=0, column=0, sticky="ew")
        self.status = tk.StringVar(value="not started")
        ttk.Label(st, textvariable=self.status, font=mono, justify="left").pack(anchor="w")
        mf = ttk.LabelFrame(side, text="Latest decoded MIB", padding=4)
        mf.grid(row=1, column=0, sticky="ew", pady=(6, 0))
        self.mib_text = tk.StringVar(value="-")
        ttk.Label(mf, textvariable=self.mib_text, font=mono, justify="left").pack(anchor="w")
        cf = ttk.LabelFrame(side, text="Cells", padding=4)
        cf.grid(row=2, column=0, sticky="ew", pady=(6, 0))
        cols = ("show", "pci", "ssbs", "mibs", "sfn", "snr", "cfo")
        self.cells = ttk.Treeview(cf, columns=cols, show="headings", height=4, selectmode="none")
        for c, w, t in zip(cols, (40, 45, 55, 55, 50, 65, 80), ("Show", "PCI", "SSBs", "MIBs", "SFN", "SNR dB", "CFO kHz")):
            self.cells.heading(c, text=t)
            self.cells.column(c, width=w, anchor="center" if c == "show" else "e")
        self.cells.pack(fill="x")
        self.cells.bind("<ButtonRelease-1>", self.on_cell_click)
        self.cells.bind("<Double-Button-1>", self.on_cell_double)
        row = ttk.Frame(cf)
        row.pack(fill="x", pady=(3, 0))
        ttk.Button(row, text="Show all", command=self.show_all).pack(side="left")
        ttk.Label(row, text="click: show/hide   double-click: only this", foreground="#555").pack(side="left", padx=8)
        lf = ttk.LabelFrame(side, text="Log", padding=4)
        lf.grid(row=3, column=0, sticky="nsew", pady=(6, 0))
        self.log = tk.Text(lf, height=10, width=44, font=("DejaVu Sans Mono", 8), wrap="none")
        sb = ttk.Scrollbar(lf, command=self.log.yview)
        self.log.configure(yscrollcommand=sb.set)
        self.log.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")

        root.protocol("WM_DELETE_WINDOW", self.close)
        root.after(100, self.poll)
        root.after(REDRAW_MS, self.redraw_tick)
        if args.autostart:
            root.after(200, self.start)
        if args.quit_after:
            root.after(int(args.quit_after * 1000), self.snapshot_and_quit)

    # ------------------------------------------------------------------ state
    def reset_data(self):
        self.cfg: dict = {}
        self.t_ref = None
        self.n_frames = self.n_iq = self.n_skipped = self.n_decoded = self.n_mib = self.n_reconnects = 0
        self.modes = {"A": 0, "T": 0}
        self.recent = deque(maxlen=200)              # host times of frames (rate)
        self.decode_ms = deque(maxlen=50)
        self.series: dict[int, dict] = {}
        self.chip = {"t": deque(maxlen=5000), "cfo": deque(maxlen=5000)}
        self.last: dict[int, deque] = {}             # latest decoded results per PCI (grid etc.)
        self.mib_text_by_pci: dict[int, tuple] = {}  # (arrival, text) of the latest MIB per PCI
        self.arrival = 0
        self.grid_pow: dict[int, deque] = {}         # |Y|^2 / noise of the last 25 SSBs per PCI
        self.colors: dict[int, str] = {}
        # the board's clock restarts at 0 when the board resets: times are continued with an
        # offset, and the SFN timing fit starts again (epoch)
        self.epoch, self.epoch_offset_us = 0, 0.0
        self.last_chip_us = self.last_host = None

    def describe_args(self) -> str:
        a = self.args
        if a.replay:
            return f"replay {Path(a.replay).name}"
        return f"live {a.port}  LO {a.freq} MHz  SSB {a.ssb} MHz  {a.rate:g} MS/s"

    def color(self, pci: int) -> str:
        if pci not in self.colors:
            self.colors[pci] = PCI_COLORS[len(self.colors) % len(PCI_COLORS)]
        return self.colors[pci]

    def t_rel(self, t_us):
        if t_us is None:
            return None
        if self.t_ref is None:
            self.t_ref = t_us
        return (t_us - self.t_ref) / 1e6

    # ------------------------------------------------------------------ control
    def start(self):
        if self.source and self.source.is_alive():
            return
        self.reset_data()
        self.log.delete("1.0", "end")
        busy = lambda: self.inflight + self.events.qsize() >= MAX_BACKLOG  # noqa: E731
        self.source = (ReplaySource(self.args, self.events, self.pool, busy) if self.args.replay
                       else LiveSource(self.args, self.events, self.pool))
        self.source.start()
        self.btn_start.configure(state="disabled")
        self.btn_stop.configure(state="normal")
        self.t_start = time.monotonic()

    def stop(self):
        if self.source:
            self.source.stop_flag.set()
        self.btn_stop.configure(state="disabled")

    def close(self):
        self.stop()
        if self.source:
            self.source.join(timeout=8)
        self.pool.shutdown(wait=False, cancel_futures=True)
        self.root.destroy()

    def save_png(self, path: Path | None = None):
        path = path or ROOT / "results" / f"ssb_live_{time.strftime('%Y%m%d_%H%M%S')}.png"
        self.fig.savefig(path, dpi=150)
        self.add_log(f"saved {path.relative_to(ROOT) if path.is_relative_to(ROOT) else path}")

    def snapshot_and_quit(self):
        if self.args.snapshot:
            self.redraw()
            self.root.update()
            self.fig.savefig(self.args.snapshot, dpi=110)
            try:
                from PIL import ImageGrab
                x, y = self.root.winfo_rootx(), self.root.winfo_rooty()
                w, h = self.root.winfo_width(), self.root.winfo_height()
                ImageGrab.grab(bbox=(x, y, x + w, y + h)).save(str(Path(self.args.snapshot).with_suffix("")) + "_window.png")
            except Exception as e:  # no X screen grab available
                print(f"window grab failed: {e}")
        print(f"ssb_live: {self.n_frames} frames, {self.n_decoded} decoded, {self.n_mib} MIBs, "
              f"{self.n_skipped} skipped" + "".join(f"; PCI {p}: {len(s['t'])} SSBs, {sum(s['ok'])} MIBs"
                                                    for p, s in self.series.items()))
        self.close()

    # ------------------------------------------------------------------ event handling
    def add_log(self, text: str):
        self.log.insert("end", text + "\n")
        if int(self.log.index("end-1c").split(".")[0]) > 400:
            self.log.delete("1.0", "50.0")
        self.log.see("end")

    def poll(self):
        try:
            while True:
                ev = self.events.get_nowait()
                k = ev["kind"]
                if k == "frame":
                    self.on_frame(ev)
                elif k == "config":
                    self.cfg = ev
                    self.src_text.set(f"{ev['source']}  LO {ev['f_lo'] / 1e6:.0f} MHz  SSB {ev['f_ssb'] / 1e6:.2f} MHz  "
                                      f"{ev['rate']:g} MS/s  gain {ev['gain']}  N_ID2 {ev['nid2']}"
                                      + (f"  recording -> {ev['record']}" if ev.get("record") else ""))
                elif k == "reconnect":
                    self.n_reconnects += 1
                elif k == "log":
                    self.add_log(ev["text"])
                elif k == "error":
                    self.add_log("ERROR " + ev["text"])
                elif k == "stopped":
                    self.btn_start.configure(state="normal")
                    self.btn_stop.configure(state="disabled")
        except queue.Empty:
            pass
        try:
            while True:
                item = self.results.get_nowait()
                if isinstance(item, Exception):
                    self.add_log(f"decode error: {type(item).__name__}: {item}")
                    continue
                for r in item:
                    self.on_result(r)
        except queue.Empty:
            pass
        self.update_status()
        self.root.after(50, self.poll)

    def continue_clock(self, ev):
        tc = ev.get("t_pss_chip_us") if ev.get("t_pss_chip_us") is not None else ev.get("t_first_us")
        if tc is None:
            return
        if self.last_chip_us is not None and tc < self.last_chip_us - 1e6:
            self.epoch += 1
            self.epoch_offset_us += self.last_chip_us - tc + max(0.0, ev["host_time"] - self.last_host) * 1e6
            for s in self.series.values():
                for k in ("tf", "sfn", "tfr"):
                    s[k].clear()
            self.add_log("board clock restarted (reset): time axis continued, timing fit restarted")
        self.last_chip_us, self.last_host = tc, ev["host_time"]
        for k in ("t_pss_chip_us", "t_first_us"):
            if ev.get(k) is not None:
                ev[k] += self.epoch_offset_us

    def on_frame(self, ev):
        self.continue_clock(ev)
        self.n_frames += 1
        self.recent.append(time.monotonic())
        self.modes[ev["mode"]] = self.modes.get(ev["mode"], 0) + 1
        t = self.t_rel(ev.get("t_pss_chip_us") if ev.get("t_pss_chip_us") is not None else ev.get("t_first_us"))
        if ev.get("cfo") is not None and t is not None and ev.get("t_pss") is not None:
            self.chip["t"].append(t)
            self.chip["cfo"].append(ev["cfo"] / 1e3)
        if ev.get("iq") is None:
            return
        self.n_iq += 1
        if not self.args.replay and self.inflight >= MAX_BACKLOG:
            self.n_skipped += 1
            return
        job = {k: ev.get(k) for k in ("iq", "fs", "t_pss", "nid2", "cfo", "f_ssb", "f_lo", "t_first_us", "seq", "mode", "stat",
                                      "R", "cfo_ddc", "pci_chip", "scale")}
        job["epoch"] = self.epoch
        with self.lock:
            self.inflight += 1
        self.pool.submit(decode_job, job).add_done_callback(self.on_done)

    def on_done(self, fut):
        """Pool thread: hand the result to the window, free the slot at once."""
        try:
            self.results.put(fut.result())
        except Exception as e:  # includes cancellation at shutdown
            self.results.put(e)
        with self.lock:
            self.inflight -= 1

    def on_result(self, r):
        self.decode_ms.append(r["decode_ms"])
        if r.get("pci") is None:
            return
        self.n_decoded += 1
        pci = r["pci"]
        s = self.series.setdefault(pci, {k: deque(maxlen=5000) for k in ("t", "cfo", "snr", "ok", "tf", "sfn", "tfr")})
        t = self.t_rel(r.get("t_pss_us"))
        if t is None:
            t = time.monotonic() - self.t_start
        s["t"].append(t)
        s["cfo"].append(r["cfo"] / 1e3)
        s["snr"].append(r["snr_db"])
        s["ok"].append(bool(r["ok"]))
        self.arrival += 1
        r["arrival"] = self.arrival
        self.last.setdefault(pci, deque(maxlen=6)).append(r)
        if r.get("grid") is not None and r.get("nvar"):
            self.grid_pow.setdefault(pci, deque(maxlen=25)).append(np.abs(r["grid"].T) ** 2 / r["nvar"])
        if r["ok"]:
            self.n_mib += 1
            m = r["mib"]
            if r.get("t_frame_us") is not None and r.get("epoch", 0) == self.epoch:
                s["tf"].append(r["t_frame_us"])
                s["sfn"].append(m["sfn"])
                s["tfr"].append(t)
            self.mib_text_by_pci[pci] = (self.arrival,
                f"PCI {pci}   SSB index {r['issb']}   half frame {m['hrf']}\n"
                f"SFN {m['sfn']:4d}         t = {t:8.3f} s\n"
                f"subCarrierSpacingCommon  {m['scs_common_khz']} kHz\n"
                f"ssb-SubcarrierOffset     k_SSB {m['kssb']}\n"
                f"dmrs-TypeA-Position      pos{m['dmrs_typeA_pos']}\n"
                f"pdcch-ConfigSIB1         {m['pdcch_config_sib1']}\n"
                f"  CORESET#0 {m['coreset0']}, searchSpace#0 {m['searchspace0']}\n"
                f"cellBarred               {'barred' if m['cell_barred'] else 'notBarred'}\n"
                f"intraFreqReselection     {m['intra_freq_reselection']}\n"
                f"SNR per RE {r['snr_db']:+5.1f} dB   CFO {r['cfo'] / 1e3:+7.3f} kHz\n"
                + (f"RSRP {r['rsrp_dbfs']:6.1f} dBFS (not calibrated)" if r.get("rsrp_dbfs") is not None else "RSRP -")
                + (f"   chip PCI {r['pci_chip']}" if r.get("pci_chip") not in (None, -1) else "") + "\n"
                f"PSS stat (chip) {r['stat'] if r['stat'] is not None else '-'}   SSS stat {r['sss_stat']:.0f}")
            if pci not in self.hidden:
                self.add_log(f"{t:8.3f}s PCI {pci} i{r['issb']} SFN {m['sfn']:4d} h{m['hrf']} "
                             f"{r['snr_db']:+5.1f} dB {r['cfo'] / 1e3:+6.2f} kHz")
        elif pci not in self.hidden:
            self.add_log(f"{t:8.3f}s PCI {pci} i{r.get('issb')} no MIB      {r['snr_db']:+5.1f} dB"
                         if r.get("snr_db") is not None else f"{t:8.3f}s PCI {pci}: {r.get('error')}")
        self.update_cell(pci)
        self.update_mib_panel()
        self.dirty = True

    # ------------------------------------------------------------------ PCI selection
    def update_cell(self, pci: int):
        d = self.series[pci]
        sfn = d["sfn"][-1] if d["sfn"] else "-"
        okc = [c for c, o in zip(d["cfo"], d["ok"]) if o]
        vals = ("✓" if pci not in self.hidden else "–", pci, len(d["t"]), sum(d["ok"]), sfn,
                f"{np.median(d['snr']):+.1f}", f"{np.median(okc):+.3f}" if okc else "-")
        iid = str(pci)
        if self.cells.exists(iid):
            self.cells.item(iid, values=vals, tags=(iid,))
        else:
            kids = [int(k) for k in self.cells.get_children()]
            self.cells.insert("", sum(1 for k in kids if k < pci), iid=iid, values=vals, tags=(iid,))
        # row text in the PCI's plot colour, grey when hidden
        self.cells.tag_configure(iid, foreground="#999999" if pci in self.hidden else self.color(pci))

    def update_mib_panel(self):
        shown = [v for p, v in self.mib_text_by_pci.items() if p not in self.hidden]
        self.mib_text.set(max(shown)[1] if shown else "-")

    def set_hidden(self, hidden: set):
        self.hidden = set(hidden)
        for p in self.series:
            self.update_cell(p)
        self.update_mib_panel()
        self.dirty = True
        if not self.paused.get():
            self.redraw()

    def on_cell_click(self, event):
        if self._skip_release:                 # second release of a double-click
            self._skip_release = False
            return
        iid = self.cells.identify_row(event.y)
        if iid:
            self.set_hidden(self.hidden ^ {int(iid)})

    def on_cell_double(self, event):
        iid = self.cells.identify_row(event.y)
        if iid:
            self._skip_release = True
            self.set_hidden(set(self.series) - {int(iid)})

    def show_all(self):
        self.set_hidden(set())

    def update_status(self):
        now = time.monotonic()
        rate = sum(1 for t in self.recent if now - t < 5) / 5
        dm = np.median(self.decode_ms) if self.decode_ms else float("nan")
        c = self.cfg
        self.status.set(
            f"frames        {self.n_frames:6d}   ({rate:4.1f} /s)\n"
            f"with I/Q      {self.n_iq:6d}   acq {self.modes.get('A', 0)} / track {self.modes.get('T', 0)}\n"
            f"decoded       {self.n_decoded:6d}   MIBs {self.n_mib}\n"
            f"skipped       {self.n_skipped:6d}   backlog {self.inflight}\n"
            f"decode time   {dm:6.1f} ms per cut   draw {self.draw_ms:4.0f} ms\n"
            f"board lost    {self.n_reconnects:6d}x (reconnected automatically)\n"
            f"SSB {c.get('f_ssb', 0) / 1e6:.2f} MHz  LO {c.get('f_lo', 0) / 1e6:.0f} MHz")

    # ------------------------------------------------------------------ plots
    def build_axes(self):
        f = self.fig
        gs = f.add_gridspec(3, 3)
        self.ax_grid = f.add_subplot(gs[0, :2])
        self.ax_con = f.add_subplot(gs[0, 2])
        self.ax_env = f.add_subplot(gs[1, :2])
        self.ax_ch = f.add_subplot(gs[1, 2])
        self.ax_cfo = f.add_subplot(gs[2, 0])
        self.ax_tim = f.add_subplot(gs[2, 1])
        self.ax_snr = f.add_subplot(gs[2, 2])
        self.im = self.ax_grid.imshow(np.full((4, 240), -10.0), aspect="auto", origin="upper", cmap="viridis",
                                      vmin=-10, vmax=20, extent=(-0.5, 239.5, 3.5, -0.5), interpolation="nearest")
        cb = f.colorbar(self.im, ax=self.ax_grid, pad=0.01)
        cb.set_label("dB over noise", fontsize=8)
        self.ax_grid.set_yticks(range(4), ["0  PSS", "1  PBCH", "2  SSS + PBCH", "3  PBCH"], fontsize=8)
        self.ax_grid.set_xlabel("subcarrier k (SS_REF = 120)", fontsize=8)
        self.ax_grid.set_title("SSB resource grid |Y|² / noise", fontsize=9, loc="left")
        q = np.array([1 + 1j, 1 - 1j, -1 + 1j, -1 - 1j]) / np.sqrt(2)    # order of sent_idx
        self.ax_con.axhline(0, color="#aaa", lw=0.7)
        self.ax_con.axvline(0, color="#aaa", lw=0.7)
        self.con = self.ax_con.scatter([], [], s=5, alpha=0.4, linewidths=0)
        self.con_mean = self.ax_con.scatter(q.real, q.imag, s=70, c=SYM_COLORS, edgecolors="black", linewidths=1.0,
                                            zorder=4, alpha=0)
        self.ax_con.plot(q.real, q.imag, "+", color="black", ms=13, mew=1.3, zorder=5)
        self.ax_con.set_xlim(-2, 2)
        self.ax_con.set_ylim(-2, 2)
        self.ax_con.set_aspect("equal")
        self.ax_con.set_title("PBCH, last 6 SSBs (colour: symbol sent)", fontsize=9, loc="left")
        self.con_note = self.ax_con.text(0.02, 0.02, "", transform=self.ax_con.transAxes, fontsize=7, va="bottom",
                                         ha="left", color="#333",
                                         bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#ccc", alpha=0.85))
        spans = [(-4.69, 0, "CP"), (0, 66.67, "PSS")] + [(l * 71.35 - 4.69, l * 71.35, "") for l in (1, 2, 3)] + \
                [(71.35, 138.02, "PBCH"), (142.71, 209.38, "SSS + PBCH"), (214.06, 280.73, "PBCH")]
        for a, b, lab in spans:
            self.ax_env.axvspan(a, b, color="#ddd" if lab in ("", "CP") else "#eef3fb", zorder=0)
            if lab and lab != "CP":
                self.ax_env.text((a + b) / 2, 1.0, lab, transform=self.ax_env.get_xaxis_transform(), ha="center",
                                 va="bottom", fontsize=7, color="#555")
        (self.env_line,) = self.ax_env.plot([], [], lw=1, color=PCI_COLORS[0])
        self.ax_env.set_xlim(-30, 320)
        self.ax_env.set_xlabel("time from PSS useful start (µs)", fontsize=8)
        self.ax_env.set_ylabel("dB (counts²)", fontsize=8)
        self.ax_env.set_title("Latest cut: power, 0.5 µs averages", fontsize=9, loc="left", pad=12)
        self.ch_lines = [self.ax_ch.plot([], [], lw=1, color="#bbb")[0] for _ in range(5)]
        (self.ch_line,) = self.ax_ch.plot([], [], lw=1.5, color=PCI_COLORS[0])
        self.ax_ch.set_xlim(0, 239)
        self.ax_ch.set_xlabel("subcarrier k", fontsize=8)
        self.ax_ch.set_ylabel("|H| (dB)", fontsize=8)
        self.ax_ch.set_title("Channel estimate (grey: previous)", fontsize=9, loc="left")
        for ax, yl, t in ((self.ax_cfo, "CFO (kHz)", "CFO"), (self.ax_tim, "residual (µs)", "Frame timing (SFN-anchored)"),
                          (self.ax_snr, "SNR per RE (dB)", "SNR (o MIB ok, x failed)")):
            ax.set_xlabel("time (s)", fontsize=8)
            ax.set_ylabel(yl, fontsize=8)
            ax.set_title(t, fontsize=9, loc="left")
        (self.chip_line,) = self.ax_cfo.plot([], [], lw=0.8, color="#999", label="chip")
        self.ax_cfo.set_title("CFO of MIB-decoded SSBs", fontsize=9, loc="left")
        self.ax_tim.set_title("Frame timing (SFN-anchored)", fontsize=9, loc="left")
        kw = dict(fontsize=7.5, va="top", ha="left", color="#333",
                  bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#ccc", alpha=0.85))
        self.cfo_note = self.ax_cfo.text(0.02, 0.98, "", transform=self.ax_cfo.transAxes, **kw)
        self.tim_note = self.ax_tim.text(0.02, 0.98, "", transform=self.ax_tim.transAxes, **kw)
        self.cfo_pts, self.tim_pts, self.snr_ok, self.snr_bad = {}, {}, {}, {}
        for ax in self.fig.axes:
            ax.tick_params(labelsize=7)
            ax.grid(True, color="#eee", lw=0.6)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)

    def redraw_tick(self):
        if self.dirty and not self.paused.get():
            self.redraw()
            self.dirty = False
        self.root.after(REDRAW_MS, self.redraw_tick)

    def redraw(self):
        t0 = time.perf_counter()
        self._redraw()
        if self.fig.get_layout_engine() is not None and self.n_decoded > 0:
            self.canvas.draw()
            self.fig.set_layout_engine("none")       # keep the positions: constrained layout costs every draw
        self.draw_ms = 1e3 * (time.perf_counter() - t0)

    def visible_last(self, n: int = 6) -> list:
        """The latest n decoded results of the PCIs that are shown, in arrival order."""
        res = [x for p, d in self.last.items() if p not in self.hidden for x in d]
        return sorted(res, key=lambda x: x["arrival"])[-n:]

    def _redraw(self):
        last = self.visible_last()
        if last:
            r = last[-1]
            col = self.color(r["pci"])
            gp = self.grid_pow.get(r["pci"])
            if gp:
                # one SSB at -5 dB per RE is mostly noise: the mean of the last ones shows the structure
                self.im.set_data(10 * np.log10(np.mean(gp, axis=0) + 1e-6))
                self.ax_grid.set_title(f"SSB resource grid |Y|² / noise, mean of the last {len(gp)} SSBs of PCI "
                                       f"{r['pci']} (latest: SSB {r.get('issb')}"
                                       + (f", SFN {r['mib']['sfn']})" if r.get("mib") else ", no MIB)"), fontsize=9, loc="left")
            done = [x for x in last if x.get("eq") is not None]
            if done:
                # MMSE-equalised symbols; after a decoded MIB each point is coloured by the symbol that
                # was sent (re-encoded), grey if the MIB failed. At -2 dB per RE the four clouds overlap
                # (noise larger than signal), but their means sit on the four QPSK points.
                e = np.concatenate([x["eq"] for x in done])
                cols = np.concatenate([np.array(SYM_COLORS)[x["sent_idx"]] if x.get("sent_idx") is not None
                                       else np.full(len(x["eq"]), "#999999") for x in done])
                self.con.set_offsets(np.c_[e.real, e.imag])
                self.con.set_facecolors(cols)
                known = [x for x in done if x.get("sent_idx") is not None]
                if known:
                    ek = np.concatenate([x["eq"] for x in known])
                    ik = np.concatenate([x["sent_idx"] for x in known])
                    means = np.array([ek[ik == k].mean() if np.any(ik == k) else np.nan for k in range(4)])
                    self.con_mean.set_offsets(np.c_[means.real, means.imag])
                    self.con_mean.set_alpha(1)
                    agree = np.mean([x["raw_agree"] for x in known])
                    self.con_note.set_text(f"raw bits right: {100 * agree:.0f} % (chance 50 %)\n"
                                           f"big dots: mean per colour; +: QPSK points\n"
                                           f"SNR per RE {np.mean([x['snr_db'] for x in known]):+.1f} dB")
            if r.get("env_t") is not None:
                self.env_line.set_data(r["env_t"], r["env_db"])
                self.env_line.set_color(col)
                self.ax_env.relim()
                self.ax_env.autoscale_view(scalex=False)
            chs = [x["channel"] for x in last if x.get("channel") is not None]
            if chs:
                k = np.arange(240)
                self.ch_line.set_data(k, 20 * np.log10(np.abs(chs[-1]) + 1e-9))
                self.ch_line.set_color(col)
                for ln, h in zip(self.ch_lines, chs[-6:-1][::-1]):
                    ln.set_data(k, 20 * np.log10(np.abs(h) + 1e-9))
                self.ax_ch.relim()
                self.ax_ch.autoscale_view(scalex=False)
        self.chip_line.set_data(list(self.chip["t"]), list(self.chip["cfo"]))
        shown = [p for p in self.series if p not in self.hidden]
        main = max(shown, key=lambda p: sum(self.series[p]["ok"]), default=None)
        if main is None:
            self.tim_note.set_text("")
            self.cfo_note.set_text("")
        for pci, s in self.series.items():
            col = self.color(pci)
            if pci not in self.cfo_pts:
                (self.cfo_pts[pci],) = self.ax_cfo.plot([], [], "o", ms=2.5, color=col, label=f"PCI {pci}")
                (self.tim_pts[pci],) = self.ax_tim.plot([], [], "o", ms=2.5, color=col, label=f"PCI {pci}")
                (self.snr_ok[pci],) = self.ax_snr.plot([], [], "o", ms=2.5, color=col, label=f"PCI {pci}")
                (self.snr_bad[pci],) = self.ax_snr.plot([], [], "x", ms=3.5, color=col)
            vis = pci not in self.hidden
            for art in (self.cfo_pts[pci], self.tim_pts[pci], self.snr_ok[pci], self.snr_bad[pci]):
                art.set_visible(vis)
            t = np.array(s["t"])
            ok = np.array(s["ok"], bool)
            # CFO of decoded SSBs only: a failed one may have slipped by the 7.0 kHz ambiguity
            self.cfo_pts[pci].set_data(t[ok], np.array(s["cfo"])[ok])
            snr = np.array(s["snr"])
            self.snr_ok[pci].set_data(t[ok], snr[ok])
            self.snr_bad[pci].set_data(t[~ok], snr[~ok])
            res = self.timing_residual(s)
            if res is not None:
                tr, rr, ppm_now, ppm_min, rms = res
                self.tim_pts[pci].set_data(tr, rr)
                if pci == main:
                    # the board's clock against the gNB's: its rate gives the CFO the crystal should cause
                    exp = -ppm_now * 1e-6 * self.cfg.get("f_ssb", 0) / 1e3
                    cf = np.array(s["cfo"])[np.array(s["ok"], bool)][-50:]
                    self.tim_note.set_text(f"PCI {pci}: {rms:.2f} µs RMS after a quadratic fit\n"
                                           f"board clock {ppm_now:+.2f} ppm now, {ppm_min:+.2f} ppm/min")
                    self.cfo_note.set_text(f"PCI {pci}, last 50: {np.median(cf):+.2f} kHz\n"
                                           f"from the clock rate: {exp:+.2f} kHz")
        # legends list the shown PCIs only
        hs = [self.cfo_pts[p] for p in sorted(shown)]
        self.ax_cfo.legend(handles=[self.chip_line] + hs, fontsize=7, loc="lower right", frameon=False)
        if len(hs) > 1:
            self.ax_snr.legend(handles=[self.snr_ok[p] for p in sorted(shown)], fontsize=7, loc="lower right", frameon=False)
        elif self.ax_snr.get_legend():
            self.ax_snr.get_legend().remove()
        for ax in (self.ax_cfo, self.ax_tim, self.ax_snr):
            ax.relim(visible_only=True)
            ax.autoscale_view()
        self.canvas.draw_idle()

    @staticmethod
    def timing_residual(s):
        """Frame start minus SFN x 10 ms (unwrapped over the 1024-frame period) on the board's
        clock, minus a quadratic fit: the board's crystal drifts as it warms up (live run of
        23 Nov: -7.19 ppm, CFO +0.5 kHz in 78 s), so a straight line leaves a +-1.5 us bow.
        Returns times, residuals (us), clock rate at the last point (ppm), its change (ppm/min), RMS."""
        if len(s["tf"]) < 4:
            return None
        tf = np.array(s["tf"])
        t0 = tf - np.array(s["sfn"]) * 1e4
        t0 = t0 + np.round((t0[0] - t0) / 10.24e6) * 10.24e6
        x = tf - tf[0]
        deg = 2 if x[-1] > 5e6 else 1
        p = np.polyfit(x, t0, deg)
        r = t0 - np.polyval(p, x)
        d = np.polyder(p)
        ppm_now = np.polyval(d, x[-1]) * 1e6
        ppm_min = (2 * p[0] * 60e6 * 1e6) if deg == 2 else 0.0
        return np.array(s["tfr"]), r, ppm_now, ppm_min, float(np.std(r))


# ============================================================================ main
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


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_source_args(p)
    p.add_argument("--hide", type=int, nargs="*", default=[], help="PCIs to hide at the start (click in the Cells table to change)")
    p.add_argument("--no-autostart", dest="autostart", action="store_false")
    p.add_argument("--geometry", default="1560x920")
    p.add_argument("--dpi", type=int, default=90)
    p.add_argument("--quit-after", type=float, default=None, help="seconds; then save --snapshot and quit (tests)")
    p.add_argument("--snapshot", default=None, help="PNG path for --quit-after")
    a = p.parse_args(argv)
    check_source_args(p, a)
    import signal
    import tkinter as tk
    root = tk.Tk()
    app = App(root, a)
    signal.signal(signal.SIGTERM, lambda *_: root.after(0, app.close))   # runs at the next poll (50 ms)
    root.mainloop()


if __name__ == "__main__":
    main()
