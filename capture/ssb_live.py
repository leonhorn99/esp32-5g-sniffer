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

Channel panel (live): type an LO and an SSB frequency (or 'auto') and press Track; or Scan a list
of LOs ('2131 2670', '2175-2715/12', or '2g 5g' for the C5's whole 2.13-2.73 and 4.79-5.99 GHz
in 12 MHz steps): full bursts at each LO, timed until every SSB position of a 20 ms period was
seen 'looks per SSB' times, are searched over the whole SSB raster on the PC
(capture/cell_scan.py: 15 kHz SSBs below 3 GHz, 30 kHz near TDD bands and at 5 GHz), and a
double-click on a found 15 kHz cell tracks it (LO, SSB, N_ID2 and CFO from the scan). The window
shows the scan's estimated duration (both ranges with 2 looks: 14 min on the board). Once
the SFN timing or a decoded cell of the scan has measured the board's crystal, scans and new
channels centre their CFO search on it (+-2 ppm instead of +-12 ppm).

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
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from .nr_mib import content_key
from .ssb_live_channel import ChannelPanel
from .ssb_live_panels import CellPanels
from .ssb_live_plots import PCI_COLORS, REDRAW_MS, Plots
from .ssb_sources import (MAX_BACKLOG, ROOT, SCAN_BURSTS, LiveSource, ReplaySource, ScanSource, Source, _worker_init,
                          add_source_args, check_source_args, decode_job)


# ============================================================================ the window
class App(CellPanels, ChannelPanel, Plots):
    """The window: state, control (Start/Stop, sources) and event handling. The panels and plots
    are mixins: ssb_live_panels (Cells, MIB), ssb_live_channel (Channel tab), ssb_live_plots."""

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
        self.run_id = 0                    # decode results of an earlier run (other channel) are dropped
        self.clock_ppm = None              # the board's crystal from the last SFN timing fit (ppm)
        self.found_cells: dict[str, tuple] = {}   # scan table row -> (LO, cell)
        self.scan_bursts = SCAN_BURSTS
        self.reset_data()

        root.title("SSB live view (ESP32-C5)")
        if args.geometry:
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

        # side panel: as wide as its contents need at the desktop's font size (a fixed width cut it
        # off with larger fonts), in tabs so that it fits the height of a 1080-line screen
        side = ttk.Frame(root, padding=(6, 4))
        side.grid(row=1, column=1, sticky="ns")
        side.rowconfigure(0, weight=1)
        self.tabs = ttk.Notebook(side)
        self.tabs.grid(row=0, column=0, sticky="nsew")
        tab_cell, tab_mib, tab_ch, tab_log = (ttk.Frame(self.tabs, padding=4) for _ in range(4))
        self.tabs.add(tab_cell, text="Cell")
        self.tabs.add(tab_mib, text="MIB")
        self.tabs.add(tab_ch, text="Channel")
        self.tabs.add(tab_log, text="Log")
        self.tab_cell, self.tab_mib, self.tab_ch, self.tab_log = tab_cell, tab_mib, tab_ch, tab_log
        self.ticker = tk.StringVar(value="")          # the latest message, below the tabs
        ttk.Label(side, textvariable=self.ticker, foreground="#555", width=1).grid(row=1, column=0, sticky="ew",
                                                                                  pady=(4, 0))
        mono = ("DejaVu Sans Mono", 9)
        st = ttk.LabelFrame(tab_cell, text="Stream", padding=4)
        st.grid(row=0, column=0, sticky="ew")
        self.status = tk.StringVar(value="not started")
        ttk.Label(st, textvariable=self.status, font=mono, justify="left").pack(anchor="w")
        mf = ttk.LabelFrame(tab_cell, text="Latest decoded MIB", padding=4)
        mf.grid(row=1, column=0, sticky="ew", pady=(6, 0))
        self.mib_text = tk.StringVar(value="-")
        ttk.Label(mf, textvariable=self.mib_text, font=mono, justify="left").pack(anchor="w")
        cf = ttk.LabelFrame(tab_cell, text="Cells", padding=4)
        cf.grid(row=2, column=0, sticky="ew", pady=(6, 0))
        cols = ("show", "pci", "ssbs", "mibs", "sfn", "snr", "cfo")
        self.cells = ttk.Treeview(cf, columns=cols, show="headings", height=4, selectmode="none")
        self.fit_columns(self.cells, ("Show", "PCI", "SSBs", "MIBs", "SFN", "SNR dB", "CFO kHz"),
                         ("✓", "1007", "99999", "99999", "1023", "-10.0", "+15.000"), center=("show",))
        self.cells.pack(fill="x")
        self.cells.bind("<ButtonRelease-1>", self.on_cell_click)
        self.cells.bind("<Double-Button-1>", self.on_cell_double)
        self.cells.bind("<Button-3>", self.on_cell_menu)
        row = ttk.Frame(cf)
        row.pack(fill="x", pady=(3, 0))
        ttk.Button(row, text="Show all", command=self.show_all).pack(side="left")
        ttk.Label(row, text="click: show/hide  double: only this\nright-click: track this PCI on the chip",
                  foreground="#555", justify="left").pack(side="left", padx=8)
        self.build_mib_panel(tab_mib)
        self.build_channel_panel(tab_ch)
        self.log = tk.Text(tab_log, height=10, width=44, font=("DejaVu Sans Mono", 8), wrap="none")
        sb = ttk.Scrollbar(tab_log, command=self.log.yview)
        self.log.configure(yscrollcommand=sb.set)
        self.log.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        if not args.geometry:
            # wide enough for the side panel next to a readable figure, tall enough for the panel
            root.update_idletasks()
            w = max(1560, side.winfo_reqwidth() + 1000)
            h = max(920, side.winfo_reqheight() + bar.winfo_reqheight() + 16)
            root.geometry(f"{w}x{h}")

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
        self.mib_info: dict[int, dict] = {}          # per PCI: MIB count, content changes, SSB indices (MIB tab)
        self.mib_dirty = True
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
    def drop_stale_events(self):
        """Before a new stream or scan: the ended source's queued frames, logs and its 'stopped'
        belong to the old channel. Scan results are kept for the table."""
        keep = []
        try:
            while True:
                ev = self.events.get_nowait()
                if ev["kind"] == "scan":
                    keep.append(ev)
        except queue.Empty:
            pass
        for ev in keep:
            self.events.put(ev)

    def start(self):
        if self.source and self.source.is_alive():
            return
        self.drop_stale_events()
        self.run_id += 1
        self.reset_data()
        # the previous run's PCIs (another channel, perhaps): table rows, plotted points, MIB panel
        for iid in self.cells.get_children():
            self.cells.delete(iid)
        for d in (self.cfo_pts, self.tim_pts, self.snr_ok, self.snr_bad):
            for art in d.values():
                art.remove()
            d.clear()
        self.mib_text.set("-")
        self.mib_dirty = True
        self.dirty = True
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
    def add_log(self, text: str, ticker: bool = True):
        if ticker:
            self.ticker.set(text if len(text) <= 90 else text[:87] + "...")
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
                elif k == "scan":
                    self.on_scan(ev)
                elif k == "stopped":
                    self.btn_start.configure(state="normal")
                    self.btn_stop.configure(state="disabled")
                    if isinstance(self.source, ScanSource):
                        self.add_log(f"scan finished: {len(self.found_cells)} cells; double-click one to track it")
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
        job["run"] = self.run_id
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
        if r.get("run") is not None and r["run"] != self.run_id:
            return                                     # decoded after a channel switch: other cell
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
            info = self.mib_info.setdefault(pci, {"n": 0, "key": None, "changes": 0, "issb": Counter(), "mib": None,
                                                  "f_ssb": None})
            key = content_key(m)
            if info["key"] is not None and key != info["key"]:
                info["changes"] += 1
                self.add_log(f"PCI {pci}: MIB content changed (other than SFN): the cell was reconfigured")
            info["n"] += 1
            info["key"], info["mib"], info["f_ssb"] = key, m, self.cfg.get("f_ssb")
            info["issb"][r["issb"]] += 1
            self.mib_dirty = True
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
                             f"{r['snr_db']:+5.1f} dB {r['cfo'] / 1e3:+6.2f} kHz", ticker=False)
        elif pci not in self.hidden:
            self.add_log(f"{t:8.3f}s PCI {pci} i{r.get('issb')} no MIB      {r['snr_db']:+5.1f} dB"
                         if r.get("snr_db") is not None else f"{t:8.3f}s PCI {pci}: {r.get('error')}", ticker=False)
        self.update_cell(pci)
        self.update_mib_panel()
        self.dirty = True

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


# ============================================================================ main
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_source_args(p)
    p.add_argument("--hide", type=int, nargs="*", default=[], help="PCIs to hide at the start (click in the Cells table to change)")
    p.add_argument("--no-autostart", dest="autostart", action="store_false")
    p.add_argument("--geometry", default=None, help="window size WxH (default: fits the side panel, at least 1560x920)")
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
