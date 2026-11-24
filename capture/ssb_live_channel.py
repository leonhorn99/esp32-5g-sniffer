"""Live view: the Channel tab (LO / SSB entry with Track, a scan over LOs, the table of cells
found). A mixin of ssb_live.App, which holds the state and starts the sources.
"""
from __future__ import annotations

from .ssb_sources import (C5_TUNE_MHZ, SCAN_LOOKS, SSB_HALF_SPAN_MHZ, ScanSource, lo_for_ssb, lo_tunable,
                          nr_band, parse_los)


class ChannelPanel:
    """Mixin of ssb_live.App: Channel tab, Track and Scan."""

    # ------------------------------------------------------------------ channel selection
    def build_channel_panel(self, parent):
        """LO / SSB entry with Track, a scan over LOs, and the table of cells found (double-click:
        track that cell). Live mode only: the board has one receiver, so a scan or a new channel
        first stops the running stream."""
        tk, ttk = self.tk, self.ttk
        a = self.args
        ch = ttk.LabelFrame(parent, text="Channel", padding=4)
        ch.pack(fill="both", expand=True)
        r1 = ttk.Frame(ch)
        r1.pack(fill="x")
        self.v_lo = tk.StringVar(value=str(a.freq or ""))
        self.v_ssb = tk.StringVar(value="auto" if a.ssb in (None, "auto") else f"{float(a.ssb):g}")
        ttk.Label(r1, text="LO").pack(side="left")
        e_lo = ttk.Entry(r1, textvariable=self.v_lo, width=6)
        e_lo.pack(side="left", padx=(3, 2))
        ttk.Label(r1, text="MHz  SSB").pack(side="left")
        e_ssb = ttk.Entry(r1, textvariable=self.v_ssb, width=9)
        e_ssb.pack(side="left", padx=(3, 2))
        ttk.Label(r1, text="MHz").pack(side="left")
        self.btn_track = ttk.Button(r1, text="Track", command=self.track_entry, width=7)
        self.btn_track.pack(side="right")
        r2 = ttk.Frame(ch)
        r2.pack(fill="x", pady=(3, 0))
        ttk.Label(r2, text="Scan LOs").pack(side="left")
        self.v_scan = tk.StringVar(value=str(a.freq or "2131"))
        e_scan = ttk.Combobox(r2, textvariable=self.v_scan, width=22,
                              values=("2g 5g", "2g", "5g", "2131 2670", str(a.freq or "2131")))
        e_scan.pack(side="left", padx=(3, 2))
        self.btn_scan = ttk.Button(r2, text="Scan", command=self.scan_entry, width=7)
        self.btn_scan.pack(side="right")
        r3 = ttk.Frame(ch)
        r3.pack(fill="x", pady=(3, 0))
        ttk.Label(r3, text="looks per SSB").pack(side="left")
        self.v_looks = tk.StringVar(value=str(SCAN_LOOKS))
        e_looks = ttk.Spinbox(r3, from_=1, to=4, textvariable=self.v_looks, width=3)
        e_looks.pack(side="left", padx=(3, 8))
        self.v_eta = tk.StringVar(value="")
        ttk.Label(r3, textvariable=self.v_eta, foreground="#555").pack(side="left")
        self.v_scan.trace_add("write", lambda *_: self.update_eta())
        self.v_looks.trace_add("write", lambda *_: self.update_eta())
        self.update_eta()
        cols = ("band", "lo", "ssb", "pci", "mibs", "snr", "cfo")
        self.found = ttk.Treeview(ch, columns=cols, show="headings", height=10, selectmode="browse")
        self.fit_columns(self.found, ("Band", "LO", "SSB MHz", "PCI", "MIBs", "SNR dB", "CFO kHz"),
                         ("n41/n38 30k", "5990", "5989.990", "1007", "99/99", "-10.0", "+15.00"))
        self.found.pack(fill="both", expand=True, pady=(3, 0))
        self.found.bind("<Double-Button-1>", self.on_found_double)
        hint = ("double-click a cell to track it\nLOs: 2131 2670, a range 2175-2715/12, or 2g / 5g\n"
                "(the C5's whole 2.13-2.73 / 4.79-5.99 GHz)" if not a.replay else "replay: no board to tune")
        ttk.Label(ch, text=hint, foreground="#555", justify="left").pack(anchor="w", pady=(2, 0))
        e_lo.bind("<Return>", lambda _e: self.track_entry())
        e_ssb.bind("<Return>", lambda _e: self.track_entry())
        e_scan.bind("<Return>", lambda _e: self.scan_entry())
        if a.replay:
            for w in (e_lo, e_ssb, e_scan, e_looks, self.btn_track, self.btn_scan):
                w.configure(state="disabled")

    def scan_plan(self):
        """(tunable LOs, untunable LOs, looks) from the scan entries; ValueError if unreadable."""
        los = parse_los(self.v_scan.get())
        looks = int(self.v_looks.get())
        if not 1 <= looks <= 8:
            raise ValueError("looks per SSB: 1 to 8")
        return [lo for lo in los if lo_tunable(lo)], [lo for lo in los if not lo_tunable(lo)], looks

    def update_eta(self):
        from .cell_scan import sweep_eta_s
        try:
            los, _, looks = self.scan_plan()
        except ValueError:
            self.v_eta.set("")
            return
        t = sweep_eta_s(los, looks, getattr(self.args, "gain", "auto") == "auto", getattr(self.args, "workers", 2))
        self.v_eta.set(f"{len(los)} LOs, about {t / 60:.0f} min" if t >= 90 else f"{len(los)} LOs, about {t:.0f} s")

    def cfo_guess(self, f_mhz: float):
        """CFO (Hz) the board's crystal causes at f_mhz, from the last timing fit; None if unknown."""
        return None if self.clock_ppm is None else -self.clock_ppm * f_mhz

    def after_stop(self, fn):
        """Stop the running stream or scan, then call fn once its thread has ended (one serial port)."""
        if self.source and self.source.is_alive():
            self.source.stop_flag.set()
            self.btn_stop.configure(state="disabled")
            self.root.after(150, lambda: self.after_stop(fn))
            return
        fn()

    def track(self, lo: int, ssb, nid2: int | None = None, cfo: float | None = None, what: str = ""):
        """Retune to a channel and restart the stream. ssb: SSB centre in MHz or 'auto' (raster scan
        at this LO first). With nid2 and cfo from a scan, the chip searches only that N_ID2 and +-10 kHz."""
        a = self.args
        a.freq, a.ssb, a.nid2 = int(lo), ssb, nid2
        guess = self.cfo_guess(lo if ssb == "auto" else ssb)
        a.cfo = cfo if cfo is not None else guess
        a.track_span = 10000 if cfo is not None else None
        self.v_lo.set(str(a.freq))
        self.v_ssb.set("auto" if ssb == "auto" else f"{float(ssb):g}")
        self.src_text.set(self.describe_args())
        self.hidden = set()
        note = (f"channel: LO {a.freq} MHz, SSB {ssb if ssb == 'auto' else f'{float(ssb):.3f}'} MHz"
                + (f" ({nr_band(float(ssb))})" if ssb != "auto" else "") + (f", {what}" if what else "")
                + (f", N_ID2 {nid2}" if nid2 is not None else "")
                + (f", CFO {a.cfo / 1e3:+.2f} kHz" + (" (from the crystal)" if cfo is None else "") if a.cfo is not None else ""))

        def go():
            self.start()
            self.add_log(note)
            self.tabs.select(self.tab_cell)
        self.after_stop(go)

    def track_entry(self):
        try:
            lo = int(float(self.v_lo.get()))
            s = self.v_ssb.get().strip().lower()
            ssb = "auto" if s in ("", "auto") else float(s)
            if ssb != "auto":
                new = lo_for_ssb(ssb, lo)
                if new != lo:
                    self.add_log(f"SSB {ssb} MHz is not within 2-{SSB_HALF_SPAN_MHZ} MHz of LO {lo} MHz: using LO {new} MHz")
                    lo = new
            elif not lo_tunable(lo):
                raise ValueError(f"LO {lo} MHz: the C5 tunes {C5_TUNE_MHZ[0]} and {C5_TUNE_MHZ[1]} MHz")
        except ValueError as e:
            self.add_log(f"channel: {e}")
            return
        self.track(lo, ssb)

    def scan_entry(self):
        from .cell_scan import sweep_eta_s
        try:
            los, bad, looks = self.scan_plan()
        except ValueError:
            self.add_log("scan: give LOs in MHz, e.g. 2131 2670, 2175-2715/12 or 2g 5g, and looks per SSB 1-8")
            return
        if bad:
            self.add_log(f"scan: skipping {bad} (the C5 tunes {C5_TUNE_MHZ[0]} and {C5_TUNE_MHZ[1]} MHz)")
        if not los:
            return
        for iid in self.found.get_children():
            self.found.delete(iid)
        self.found_cells.clear()
        ppm = self.clock_ppm
        eta = sweep_eta_s(los, looks, getattr(self.args, "gain", "auto") == "auto", getattr(self.args, "workers", 2))

        def go():
            self.drop_stale_events()
            self.source = ScanSource(self.args, self.events, self.pool, los, ppm, looks, self.scan_bursts)
            self.source.start()
            self.btn_start.configure(state="disabled")
            self.btn_stop.configure(state="normal")
            self.tabs.select(self.tab_ch)
            self.add_log(f"scan: {len(los)} LOs ({los[0]}-{los[-1]} MHz), each until every SSB position was seen "
                         f"{looks}x (at most {self.scan_bursts} bursts), about {eta / 60:.0f} min; "
                         + (f"CFO around the crystal ({ppm:+.2f} ppm) +-2 ppm" if ppm is not None else
                            "crystal unknown: CFO +-12 ppm until a cell decodes"))
        self.after_stop(go)

    def on_scan(self, ev):
        lo = ev["lo"]
        if ev.get("ppm") is not None:
            self.clock_ppm = float(ev["ppm"])                  # the crystal from the scan's decoded cells
        if not ev["cells"]:
            self.add_log(f"scan LO {lo} MHz: no cell", ticker=False)
        for c in ev["cells"]:
            scs = c.get("scs_khz", 15)
            iid = f"{lo}_{c['f_ssb_mhz']:.3f}_{c['pci']}" + ("" if scs == 15 else f"_{scs:g}")
            cfo = f"{c['cfo_hz'] / 1e3:+.2f}" if c.get("cfo_hz") is not None else "-"
            band = nr_band(c["f_ssb_mhz"]) + ("" if scs == 15 else f" {scs:g}k")
            vals = (band, lo, f"{c['f_ssb_mhz']:.3f}", c["pci"], f"{c['mibs']}/{c['ssbs']}", f"{c['snr_db']:+.1f}", cfo)
            if self.found.exists(iid):
                self.found.item(iid, values=vals)
            else:
                self.found.insert("", "end", iid=iid, values=vals)
            self.found_cells[iid] = (lo, c)
            self.add_log(f"scan LO {lo} MHz: PCI {c['pci']} at {c['f_ssb_mhz']:.3f} MHz ({band}), "
                         f"{c['mibs']} MIBs of {c['ssbs']} SSBs, SNR {c['snr_db']:+.1f} dB")
        if "n" in ev:                                          # progress of a scan, in the status line
            left = ev.get("left_s") or 0.0
            self.ticker.set(f"scan: LO {ev['index'] + 1}/{ev['n']} ({lo} MHz), {len(self.found_cells)} cells, "
                            f"about {left / 60:.0f} min left")

    def on_found_double(self, event):
        iid = self.found.identify_row(event.y)
        if iid not in self.found_cells:
            return
        lo, c = self.found_cells[iid]
        if c.get("scs_khz", 15) != 15:
            self.add_log(f"PCI {c['pci']} sends 30 kHz SSBs: the chip's detector and the live decoder handle 15 kHz only")
            return
        self.track(lo, round(c["f_ssb_mhz"], 3), nid2=c["pci"] % 3, cfo=c.get("cfo_hz"), what=f"PCI {c['pci']}")
