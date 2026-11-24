"""Live view: the Cells table (show, hide, track a PCI) and the MIB tab. A mixin of ssb_live.App,
which holds the state (self.series, self.hidden, ...) and the Tk widgets.
"""
from __future__ import annotations

import numpy as np

from .nr_mib import FIELDS as MIB_FIELDS, HELP as MIB_HELP, describe as describe_mib


class CellPanels:
    """Mixin of ssb_live.App: Cells table and MIB tab."""

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

    MIB_COLS = 3                                     # cells side by side in the MIB tab

    def build_mib_panel(self, parent):
        """The MIB of up to three cells side by side, every field in words (capture/nr_mib.py): click a
        row for what it means. Below, what follows from it (sectors of one site, reconfigurations)."""
        tk, ttk = self.tk, self.ttk
        f = ttk.LabelFrame(parent, text="MIB per cell (read from the PBCH)", padding=4)
        f.pack(fill="both", expand=True)
        cols = ("field",) + tuple(f"c{j}" for j in range(self.MIB_COLS))
        self.mib_tree = ttk.Treeview(f, columns=cols, show="headings", height=len(MIB_FIELDS), selectmode="browse")
        self.fit_columns(self.mib_tree, ("Field",) + ("PCI 1007",) * self.MIB_COLS,
                         ("subCarrierSpacingCommon",) + ("2120.68-2129.32",) * self.MIB_COLS)
        self.mib_tree.column("field", anchor="w")
        for j in range(self.MIB_COLS):
            self.mib_tree.heading(f"c{j}", text="")
        for k, label, _ in MIB_FIELDS:
            self.mib_tree.insert("", "end", iid=k, values=(label,) + ("",) * self.MIB_COLS)
        self.mib_tree.pack(fill="x")
        self.mib_tree.bind("<<TreeviewSelect>>", self.on_mib_select)
        wrap = sum(int(self.mib_tree.column(c, "width")) for c in cols)
        self.mib_help = tk.StringVar(value="Click a row to see what it means.")
        ttk.Label(f, textvariable=self.mib_help, wraplength=wrap, justify="left").pack(fill="x", pady=(6, 0))
        self.mib_notes = tk.StringVar(value="")
        ttk.Label(f, textvariable=self.mib_notes, wraplength=wrap, justify="left", foreground="#555").pack(
            fill="x", pady=(6, 0))

    def on_mib_select(self, _event=None):
        sel = self.mib_tree.selection()
        if sel:
            label = next(lb for k, lb, _ in MIB_FIELDS if k == sel[0])
            self.mib_help.set(f"{label}: {MIB_HELP[sel[0]]}")

    def update_mib_tab(self):
        shown = [p for p in self.mib_info if p not in self.hidden]
        pcis = sorted(sorted(shown, key=lambda p: -self.mib_info[p]["n"])[:self.MIB_COLS])
        for j in range(self.MIB_COLS):
            col, vals, head = f"c{j}", {}, ""
            if j < len(pcis):
                pci, info = pcis[j], self.mib_info[pcis[j]]
                vals = describe_mib(info["mib"], info["f_ssb"], info["issb"].most_common(1)[0][0])
                vals.update(n=str(info["n"]), issb=", ".join(str(i) for i in sorted(info["issb"])),
                            content="constant" if not info["changes"] else f"changed {info['changes']}×")
                head = f"PCI {pci}"
            self.mib_tree.heading(col, text=head)
            for k, _, _ in MIB_FIELDS:
                self.mib_tree.set(k, col, vals.get(k, ""))
        notes, groups = [], {}
        for p in shown:
            i = self.mib_info[p]
            groups.setdefault((round((i["f_ssb"] or 0) / 1e3), i["key"]), []).append(p)
        for ps in groups.values():
            if len(ps) > 1:
                notes.append(f"PCI {' and '.join(map(str, sorted(ps)))}: same SSB frequency and the same MIB, "
                             "most likely sectors of one site.")
        ann = sorted(p for p in shown if self.mib_info[p]["mib"]["kssb"] < 24)
        if ann:
            notes.append(f"PCI {', '.join(map(str, ann))}: SIB1 announced, so phones can probably attach directly "
                         "(standalone 5G). An inference: SIB1 itself is not decoded.")
        self.mib_notes.set("\n".join(notes))

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

    def on_cell_menu(self, event):
        iid = self.cells.identify_row(event.y)
        if not iid:
            return
        pci = int(iid)
        m = self.tk.Menu(self.root, tearoff=0)
        live = not self.args.replay
        m.add_command(label=f"Track PCI {pci} on the chip (N_ID2 {pci % 3})", command=lambda: self.track_cell(pci),
                      state="normal" if live else "disabled")
        m.add_command(label=f"Show only PCI {pci}", command=lambda: self.set_hidden(set(self.series) - {pci}))
        m.tk_popup(event.x_root, event.y_root)

    def track_cell(self, pci: int):
        """The chip follows one N_ID2: switch it to this PCI's (e.g. the other sector of a site,
        seen through the PC's extra N_ID2 tries) at the same LO and SSB, from its measured CFO."""
        s = self.series.get(pci)
        okc = [c for c, o in zip(s["cfo"], s["ok"]) if o] if s else []
        cfo = float(np.median(okc)) * 1e3 if okc else None
        lo = int(round(self.cfg.get("f_lo", self.args.freq * 1e6) / 1e6))
        ssb = round(self.cfg.get("f_ssb", float(self.args.ssb) * 1e6) / 1e6, 3)
        self.track(lo, ssb, nid2=pci % 3, cfo=cfo, what=f"PCI {pci}")

    @staticmethod
    def fit_columns(tree, heads, samples, center=()):
        """Column widths from the fonts in use (the desktop may use larger ones than our tests)."""
        import tkinter.font as tkf
        body, head = tkf.nametofont("TkDefaultFont"), tkf.nametofont("TkHeadingFont")
        for c, h, smp in zip(tree["columns"], heads, samples):
            tree.heading(c, text=h)
            tree.column(c, width=max(head.measure(h) + 14, body.measure(smp) + 12), stretch=False,
                        anchor="center" if c in center else "e")
