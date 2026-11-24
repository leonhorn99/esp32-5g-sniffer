"""PCI selection in the live view (run on a virtual display by tests/test_ssb_live.py)."""
import types
import tkinter as tk

import numpy as np

from capture.ssb_live import App

args = types.SimpleNamespace(replay="x", workers=1, geometry="1560x920", dpi=90, autostart=False, quit_after=None,
                             snapshot=None, hide=[], port=None, freq=None, ssb=None, rate=20)
root = tk.Tk()
app = App(root, args)
app.t_start = 0.0


def result(pci, k):
    mib = dict(sfn=k, hrf=0, scs_common_khz=15, kssb=6, dmrs_typeA_pos=3, pdcch_config_sib1=117, coreset0=7,
               searchspace0=5, cell_barred=False, intra_freq_reselection="allowed")
    return dict(pci=pci, ok=True, cfo=15e3, snr_db=-3.0, mib=mib, issb=0, t_pss_us=1e6 + 2e4 * k,
                t_frame_us=1e6 + 2e4 * k - 148, grid=np.ones((240, 4)), nvar=1.0, eq=np.ones(432, complex),
                channel=np.ones(240), decode_ms=5.0, stat=20.0, sss_stat=30.0, epoch=0, env_t=None)


for k in range(6):
    app.on_result(result(165, k))
    app.on_result(result(166, k))
root.update()
row = lambda pci: types.SimpleNamespace(y=app.cells.bbox(str(pci))[1] + 5)  # noqa: E731
assert {x["pci"] for x in app.visible_last()} == {165, 166}
app.on_cell_click(row(165))                                   # click: hide 165
assert app.hidden == {165}, app.hidden
assert all(x["pci"] == 166 for x in app.visible_last())
assert app.mib_text.get().startswith("PCI 166")
assert app.cells.item("165")["values"][0] == "–"
app.redraw()
assert not app.cfo_pts[165].get_visible() and app.cfo_pts[166].get_visible()
app.on_cell_double(row(165))                                  # double-click: only 165
app.on_cell_click(row(165))                                   # its second release is ignored
assert app.hidden == {166}, app.hidden
assert app.mib_text.get().startswith("PCI 165")
app.show_all()
assert app.hidden == set()
app.on_cell_click(row(166))
app.on_cell_click(row(166))                                   # click twice: shown again
assert app.hidden == set()
print("PCI selection OK")
app.pool.shutdown(wait=False)
root.destroy()
