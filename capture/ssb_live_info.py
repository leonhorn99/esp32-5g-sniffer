"""Live view: the System info tab (SIB1, SIB2... in words; decode them from the window) and the
Neighbours tab (arrival times of the neighbour cells in the same cuts). A mixin of ssb_live.App, which
holds the state and the event loop.
"""
from __future__ import annotations

import threading
import time
from collections import deque

import numpy as np

from .tdoa_neighbours import Accumulator, peak

NB_BLOCK_CUTS = 250          # cuts per block for the spread of a neighbour's delay (10 s at 25 SSBs/s)
NB_DETECT = 8.0              # SSS-only peak/floor that counts as detected (capture/tdoa_neighbours.py --detect)
NB_WINDOW_US = 15.0          # +- around the serving cell (the cut holds about +-20 us)
C = 299792458.0


class InfoPanels:
    """Mixin of ssb_live.App: System info and Neighbours tabs."""

    # ------------------------------------------------------------------ state
    def info_reset(self):
        """New run (channel): neighbour sums start again; decoded SIBs stay (they are per cell)."""
        self.nb_acc: dict[int, Accumulator] = {}          # serving PCI -> sums over all cuts
        self.nb_block: dict[int, Accumulator] = {}        # serving PCI -> sums over the current block
        self.nb_taus: dict[int, dict[int, deque]] = {}    # serving PCI -> neighbour -> delay per block (us)
        self.nb_fs = None
        self.nb_dirty = True

    def init_info(self):
        self.sib_cache: dict[tuple, dict] = {}             # (pci, ssb kHz) -> nr_sib1.load_decoded()
        self.sib_shown = None
        self.sib_loading: set = set()
        self.sib_auto_done = False
        self.nb_pcis: list[int] = list(getattr(self.args, "neighbours", None) or [])
        self.nb_from = "command line" if self.nb_pcis else None
        self.nb_last_draw = 0.0

    def info_tick(self):
        """From the event loop: the Neighbours table (1/s), the decoded SIBs of the cell shown, and --sib."""
        from .ssb_sources import LiveSource, ScanSource, SibSource
        self.update_nb_tab()
        now = time.monotonic()
        if now - getattr(self, "_info_t", 0.0) < 1.0:
            return
        self._info_t = now
        self.refresh_sib()
        if getattr(self.args, "sib", False) and not self.sib_auto_done and not self.args.replay:
            tgt = self.info_target()
            src = self.source
            if tgt and self.mib_info[tgt[0]]["n"] >= 3 and isinstance(src, LiveSource) and src.is_alive() \
                    and not isinstance(src, (ScanSource, SibSource)):
                self.sib_auto_done = True
                self.add_log(f"--sib: decoding the system information of PCI {tgt[0]}")
                self.decode_sib()

    # ------------------------------------------------------------------ System info tab
    def build_sib_panel(self, parent):
        tk, ttk = self.tk, self.ttk
        f = ttk.LabelFrame(parent, text="System information (SIB1, SIB2...)", padding=4)
        f.pack(fill="both", expand=True)
        self.sib_head = tk.StringVar(value="no cell yet")
        head = ttk.Label(f, textvariable=self.sib_head, justify="left")
        head.pack(anchor="w", fill="x")
        self.sib_tree = ttk.Treeview(f, columns=("value",), show="tree headings", height=20, selectmode="browse")
        self.fit_columns(self.sib_tree, ("Value",), ("52 RB at kHz15: 2665.320-2674.680 MHz (cen",))
        self.sib_tree.heading("#0", text="Field")
        import tkinter.font as tkf
        self.sib_tree.column("#0", width=tkf.nametofont("TkDefaultFont").measure("neighbour cells, same freq. (15)xxxxxx"),
                             stretch=False)
        self.sib_tree.column("value", anchor="w")
        self.sib_tree.pack(fill="both", expand=True, pady=(3, 0))
        self.sib_tree.bind("<<TreeviewSelect>>", self.on_sib_select)
        wrap = int(self.sib_tree.column("#0", "width")) + int(self.sib_tree.column("value", "width"))
        head.configure(wraplength=wrap)
        self.sib_detail = tk.StringVar(value="Click a row for its whole text.")
        ttk.Label(f, textvariable=self.sib_detail, wraplength=wrap, justify="left").pack(fill="x", pady=(4, 0))
        row = ttk.Frame(f)
        row.pack(fill="x", pady=(4, 0))
        self.btn_sib = ttk.Button(row, text="Decode now", command=self.decode_sib,
                                  state="disabled" if self.args.replay else "normal")
        self.btn_sib.pack(side="left")
        ttk.Button(row, text="Reload", command=lambda: self.refresh_sib(force=True)).pack(side="left", padx=4)
        self.sib_status = tk.StringVar(value=("Decode now: bursts timed to the cell's SIB1 slots, then MATLAB "
                                              "(5G Toolbox) and pycrate; then the other SI messages. About 2-5 min; "
                                              "the live stream pauses meanwhile.") if not self.args.replay else
                                     "replay: decoded SIBs from results/ only")
        ttk.Label(f, textvariable=self.sib_status, wraplength=wrap, justify="left", foreground="#555").pack(
            fill="x", pady=(4, 0))

    def info_target(self):
        """(PCI, SSB MHz) of the cell the tabs describe: the shown PCI with the most MIBs."""
        shown = [p for p in self.mib_info if p not in self.hidden and self.mib_info[p].get("f_ssb")]
        if not shown:
            return None
        p = max(shown, key=lambda q: self.mib_info[q]["n"])
        return p, round(self.mib_info[p]["f_ssb"] / 1e6, 3)

    def refresh_sib(self, force: bool = False):
        """Show the decoded SIBs of the target cell; load them from results/ in a thread (pycrate is slow)."""
        tgt = self.info_target()
        if tgt is None:
            return
        key = (tgt[0], round(tgt[1] * 1e3))
        if force:
            self.sib_cache.pop(key, None)
        if key in self.sib_cache:
            if self.sib_shown != key or force:
                self.show_sib(key)
            return
        if key in self.sib_loading:
            return
        self.sib_loading.add(key)
        self.sib_head.set(f"PCI {tgt[0]} at {tgt[1]:.3f} MHz: loading decoded SIBs from results/ ...")

        def load():
            from .nr_sib1 import load_decoded
            try:
                info = load_decoded(tgt[0], tgt[1])
            except Exception as e:  # shown in the tab
                info = {"sib1": None, "sibs": {}, "files": [], "error": f"{type(e).__name__}: {e}", "plan": None}
            self.events.put({"kind": "sibinfo", "key": key, "info": info})
        threading.Thread(target=load, daemon=True).start()

    def on_sibinfo(self, ev):
        key, info = ev["key"], ev["info"]
        self.sib_loading.discard(key)
        self.sib_cache[key] = info
        if not self.nb_pcis or self.nb_from and self.nb_from.startswith("SIB3"):
            self.use_sib3(quiet=True)
        tgt = self.info_target()
        if tgt and key == (tgt[0], round(tgt[1] * 1e3)):
            self.show_sib(key)
        self.mib_dirty = True

    def show_sib(self, key):
        info = self.sib_cache[key]
        pci, ssb = key[0], key[1] / 1e3
        self.sib_shown = key
        t = self.sib_tree
        t.delete(*t.get_children())
        self.sib_rows = {}

        def add(parent, label, text):
            iid = t.insert(parent, "end", text=label, values=(text,), open=True)
            self.sib_rows[iid] = f"{label}: {text}"
            return iid
        if info.get("sib1") is not None:
            sess, n = info["sib1_from"]
            top = add("", "SIB1", f"decoded {n}x, {sess}")
            for label, text in info["sib1_rows"]:
                add(top, label, text)
            p = info.get("plan")
            if p:
                for m in p["messages"]:
                    add(top, "SI window", f"{'+'.join(m['sibs'])}: frames {m['frames']} of every {m['T']} "
                                          f"({p['window_slots']} slots), search space {p['search_space']}")
        for typ in sorted(info.get("sibs", {})):
            val, rows, sess = info["sibs"][typ]
            top = add("", typ.upper(), sess)
            for label, text in rows:
                add(top, label, text)
        n_files = len(info.get("files", []))
        if info.get("error"):
            self.sib_head.set(f"PCI {pci} at {ssb:.3f} MHz: {info['error']}")
        elif info.get("sib1") is None and not info.get("sibs"):
            self.sib_head.set(f"PCI {pci} at {ssb:.3f} MHz: nothing decoded yet ({n_files} result files). "
                              "Press Decode now (live).")
        else:
            self.sib_head.set(f"PCI {pci} at {ssb:.3f} MHz, from {n_files} result file{'s' if n_files != 1 else ''} in "
                              "results/. The TAC and cell identity locate the site: keep them out of public files.")

    def on_sib_select(self, _event=None):
        sel = self.sib_tree.selection()
        if sel and sel[0] in getattr(self, "sib_rows", {}):
            self.sib_detail.set(self.sib_rows[sel[0]])

    def sib_decoded(self, pci: int) -> bool:
        return any(k[0] == pci and v.get("sib1") is not None for k, v in self.sib_cache.items())

    def decode_sib(self):
        """Pause the stream, run SibSource on the target cell, then track it again."""
        from .sib1_capture import lo_for_coreset
        from .ssb_sources import SibSource
        if self.args.replay:
            self.add_log("SIB: a replay has no board; decode recorded bursts with matlab/sib1_offline.m")
            return
        tgt = self.info_target()
        if tgt is None:
            self.add_log("SIB: no cell with a decoded MIB yet")
            return
        pci, ssb = tgt
        info = self.mib_info[pci]
        mib = info["mib"]
        s = self.series.get(pci)
        okc = [c for c, o in zip(s["cfo"], s["ok"]) if o] if s else []
        cfo = float(np.median(okc)) * 1e3 if okc else (self.args.cfo or 0.0)
        try:
            lo = lo_for_coreset(mib, ssb * 1e6, self.args.freq)
        except ValueError as e:
            self.add_log(f"SIB: {e}")
            return
        gain = self.cfg.get("gain") if lo == self.args.freq else "auto"

        def go():
            self.drop_stale_events()
            self.source = SibSource(self.args, self.events, self.pool, pci, ssb, mib, lo, gain, cfo)
            self.source.start()
            self.btn_start.configure(state="disabled")
            self.btn_stop.configure(state="normal")
            self.btn_sib.configure(state="disabled")
            self.tabs.select(self.tab_sib)
            self.sib_resume = (pci, cfo)
            self.add_log(f"SIB: PCI {pci} at {ssb:.3f} MHz, LO {lo} MHz, CFO {cfo / 1e3:+.2f} kHz; the stream pauses")
        self.after_stop(go)

    def on_sib_event(self, ev):
        key = (ev["pci"], round(ev["ssb_mhz"] * 1e3))
        self.sib_cache.pop(key, None)
        self.refresh_sib(force=True)

    def after_sib(self):
        """SibSource ended: track the same cell again (its N_ID2 and CFO known)."""
        self.btn_sib.configure(state="normal")
        pci, cfo = getattr(self, "sib_resume", (None, None))
        if pci is not None:
            self.args.nid2, self.args.cfo, self.args.track_span = pci % 3, cfo, 10000
        self.root.after(300, self.start)

    # ------------------------------------------------------------------ Neighbours tab
    def build_nb_panel(self, parent):
        tk, ttk = self.tk, self.ttk
        f = ttk.LabelFrame(parent, text="Neighbour cells in the same cuts", padding=4)
        f.pack(fill="both", expand=True)
        r1 = ttk.Frame(f)
        r1.pack(fill="x")
        ttk.Label(r1, text="PCIs").pack(side="left")
        self.v_nb = tk.StringVar(value=" ".join(map(str, self.nb_pcis)))
        e = ttk.Entry(r1, textvariable=self.v_nb, width=30)
        e.pack(side="left", padx=(3, 3), fill="x", expand=True)
        e.bind("<Return>", lambda _e: self.apply_nb())
        r2 = ttk.Frame(f)
        r2.pack(fill="x", pady=(3, 0))
        ttk.Button(r2, text="Apply", command=self.apply_nb).pack(side="left")
        ttk.Button(r2, text="From SIB3", command=self.use_sib3).pack(side="left", padx=3)
        ttk.Button(r2, text="Reset sums", command=self.reset_nb).pack(side="left")
        self.nb_head = tk.StringVar(value="")
        ttk.Label(f, textvariable=self.nb_head, justify="left").pack(anchor="w", pady=(4, 0))
        cols = ("pci", "stat", "dt", "dm", "sd", "note")
        self.nb_tree = ttk.Treeview(f, columns=cols, show="headings", height=12, selectmode="none")
        self.fit_columns(self.nb_tree, ("PCI", "SSS stat", "Δt µs", "Δ m", "± m", ""),
                         ("1007", "999.9", "+10.000", "+3000", "999", "detected"))
        self.nb_tree.pack(fill="x", pady=(3, 0))
        wrap = sum(int(self.nb_tree.column(c, "width")) for c in cols)
        ttk.Label(f, wraplength=wrap, justify="left", foreground="#555", text=(
            "Each neighbour's SSS (and PSS) is correlated in every cut of the tracked cell, after its own PSS and "
            "SSS are subtracted, and summed. SSS stat: peak over floor, about 1-3 for noise, from 8 on detected. "
            "Δt: how much later than the tracked cell it arrives, Δ m the same as a path difference; this includes "
            "any offset between the two sites' frame timing. ± m: spread of Δ over blocks of "
            f"{NB_BLOCK_CUTS} cuts. A neighbour on another site arrives within a few µs only if the network is "
            "frame-synchronised.")).pack(fill="x", pady=(4, 0))
        self.nb_head.set(self.nb_head_text())

    def nb_head_text(self) -> str:
        if not self.nb_pcis:
            return "no neighbour list: type PCIs, or decode SIB3 (System info tab)"
        return f"{len(self.nb_pcis)} PCIs ({self.nb_from}); waiting for cuts"

    def nb_job(self):
        """What a decode job needs for the neighbours (None: no list)."""
        return {"pcis": self.nb_pcis, "window_us": NB_WINDOW_US, "cancel": True} if self.nb_pcis else None

    def apply_nb(self):
        try:
            pcis = [int(t) for t in self.v_nb.get().replace(",", " ").split()]
        except ValueError:
            self.add_log("neighbours: PCIs as numbers 0-1007, separated by spaces")
            return
        if any(not 0 <= p <= 1007 for p in pcis):
            self.add_log("neighbours: PCIs are 0-1007")
            return
        self.nb_pcis, self.nb_from = sorted(set(pcis)), "typed"
        self.reset_nb()

    def use_sib3(self, quiet: bool = False):
        from .nr_sib1 import neighbour_pcis
        tgt = self.info_target()
        info = self.sib_cache.get((tgt[0], round(tgt[1] * 1e3))) if tgt else None
        pcis = neighbour_pcis(info) if info else []
        if not pcis:
            if not quiet:
                self.add_log("neighbours: no decoded SIB3 for this cell (System info tab: Decode now)")
            return
        if pcis == self.nb_pcis:
            return
        self.nb_pcis, self.nb_from = pcis, f"SIB3 of PCI {tgt[0]}"
        self.v_nb.set(" ".join(map(str, pcis)))
        self.reset_nb()
        self.add_log(f"neighbours: {len(pcis)} PCIs from SIB3 of PCI {tgt[0]}")

    def reset_nb(self):
        self.info_reset()
        for iid in self.nb_tree.get_children():
            self.nb_tree.delete(iid)
        self.nb_head.set(self.nb_head_text())

    def on_nb(self, nb: dict):
        s = nb["serving"]
        if self.nb_fs is not None and nb["fs"] != self.nb_fs:
            self.info_reset()                       # other cut rate (raw / decimated): sums do not mix
        self.nb_fs = nb["fs"]
        curves = {p: v for p, v in nb["curves"].items() if p in self.nb_pcis}
        self.nb_acc.setdefault(s, Accumulator()).add(curves)
        blk = self.nb_block.setdefault(s, Accumulator())
        blk.add(curves)
        if blk.n >= NB_BLOCK_CUTS:
            taus = self.nb_taus.setdefault(s, {})
            for p in blk.acc:
                taus.setdefault(p, deque(maxlen=20)).append(peak(blk.acc[p] / blk.n, self.nb_fs)[1])
            self.nb_block[s] = Accumulator()
        self.nb_dirty = True

    def update_nb_tab(self):
        now = time.monotonic()
        if not self.nb_dirty or now - self.nb_last_draw < 1.0:
            return
        self.nb_last_draw, self.nb_dirty = now, False
        if not self.nb_acc:
            self.nb_head.set(self.nb_head_text())
            return
        s = max(self.nb_acc, key=lambda k: self.nb_acc[k].n)
        acc = self.nb_acc[s]
        rows = acc.rows(self.nb_fs)
        taus = self.nb_taus.get(s, {})
        n_det = sum(r["stat_sss"] >= NB_DETECT for r in rows)
        self.nb_head.set(f"relative to PCI {s}: {acc.n} cuts, {len(rows)} PCIs ({self.nb_from}), {n_det} detected")
        tree = self.nb_tree
        for iid in tree.get_children():
            tree.delete(iid)
        for r in rows:
            det = r["stat_sss"] >= NB_DETECT
            t = taus.get(r["pci"], ())
            sd = f"{np.std(t, ddof=1) * 1e-6 * C:.0f}" if det and len(t) >= 3 else ""
            note = "detected" if det else ("weak?" if r["stat_sss"] >= 5 else "")
            vals = (r["pci"], f"{r['stat_sss']:.1f}", f"{r['tau_us']:+.3f}" if det else "",
                    f"{r['delta_m']:+.0f}" if det else "", sd, note)
            tree.insert("", "end", values=vals, tags=("det",) if det else ())
        tree.tag_configure("det", foreground="#1a5fb4")
