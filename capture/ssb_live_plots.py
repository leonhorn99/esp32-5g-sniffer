"""Live view: the figure (SSB resource grid, PBCH constellation, power of the latest cut, channel
estimate, CFO, frame timing, SNR) and its periodic redraw.
A mixin of ssb_live.App, which holds the data series (self.series) and the figure (self.fig).
"""
from __future__ import annotations

import time

import numpy as np

REDRAW_MS = 700            # a redraw takes 130-300 ms of the GUI process; results queue up meanwhile
PCI_COLORS = ("#2a6fdb", "#e07b00", "#1a9e5a", "#b23a9a", "#7a6a00", "#5a5a5a")
SYM_COLORS = ("#2a6fdb", "#e07b00", "#1a9e5a", "#b23a9a")    # sent QPSK symbol (+,+) (+,-) (-,+) (-,-)


class Plots:
    """Mixin of ssb_live.App: build_axes and the periodic redraw."""

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
        if self.mib_dirty:
            self.update_mib_tab()
            self.mib_dirty = False
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
                    if len(s["tf"]) >= 20:
                        self.clock_ppm = float(ppm_now)        # centres the CFO search on other channels
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
