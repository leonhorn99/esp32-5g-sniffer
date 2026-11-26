"""System info and Neighbours tabs of the live view (run on a virtual display by tests/test_ssb_live.py).

Part 1 (always): the tabs with injected data: a decoded-SIB record (as nr_sib1.load_decoded returns it)
fills the System info tree and the neighbour list; neighbour curves from decode results fill the table.
Part 2 (private recording present): a replay of PCI 317's cuts of 26 Nov 13:30 with the neighbours of
its SIB3; PCI 100 / 90 / 69 must come out as capture/tdoa_neighbours.py found them (+1.875 / +0.567 /
+0.118 us). Optional argument: a folder for window grabs of both tabs.
"""
import sys
import time
import tkinter as tk
import types
from pathlib import Path

import numpy as np

from capture.ssb_live import App

ROOT = Path(__file__).resolve().parents[2]
SESSION = ROOT / "data" / "20261126_133006_n7_pci317_tdoa_newpos"


def make_args(**kw):
    a = dict(replay=None, port=None, freq=2670, ssb=2667.75, nid2=None, cfo=None, gain="40", rate=20, bits=8, thr=14.0,
             pre=20, post=20, period=20.0, every=1, scan_bursts=150, cfo_span=30000, raw_cuts=False, record=False,
             label="t", note="", workers=2, bandwidth=None, dccal=False, speed=0.0, hide=[], geometry="1560x980",
             dpi=90, autostart=False, quit_after=None, snapshot=None, neighbours=None, sib=False)
    a.update(kw)
    return types.SimpleNamespace(**a)


def pump(root, cond, timeout=60.0):
    t_end = time.monotonic() + timeout
    while time.monotonic() < t_end:
        root.update()
        if cond():
            return True
        time.sleep(0.02)
    return False


def grab(root, path):
    try:
        from PIL import ImageGrab
        root.update()
        x, y, w, h = root.winfo_rootx(), root.winfo_rooty(), root.winfo_width(), root.winfo_height()
        ImageGrab.grab(bbox=(x, y, x + w, y + h)).save(str(path))
    except Exception as e:  # no X screen grab available
        print(f"window grab failed: {e}")


MIB = dict(sfn=1, hrf=0, scs_common_khz=15, kssb=6, dmrs_typeA_pos=3, pdcch_config_sib1=16, coreset0=1,
           searchspace0=0, cell_barred=False, intra_freq_reselection="allowed")


def injected():
    root = tk.Tk()
    app = App(root, make_args(replay="none", ssb=2667.75))      # replay: no board, Decode now disabled
    assert str(app.btn_sib.cget("state")) == "disabled"
    app.t_start = time.monotonic()                       # as start() would set it
    app.on_result(dict(pci=317, ok=True, cfo=13.9e3, snr_db=4.0, mib=MIB, issb=3, t_pss_us=2e6, t_frame_us=2e6 - 1000,
                       grid=np.ones((240, 4)), nvar=1.0, eq=np.ones(432, complex), channel=np.ones(240), decode_ms=5.0,
                       stat=30.0, sss_stat=50.0, epoch=0, env_t=None, run=app.run_id))
    app.cfg = {"f_ssb": 2667.75e6}
    app.mib_info[317]["f_ssb"] = 2667.75e6
    assert app.info_target() == (317, 2667.75)
    sib3 = {"intraFreqNeighCellList": [{"physCellId": p} for p in (46, 69, 90, 100)]}
    info = {"sib1": {"x": 1}, "kssb": 6, "sib1_from": ("20261126_test", 3),
            "sib1_rows": [("network (PLMN) 1", "001-01"), ("bands", "n7")],
            "sibs": {"sib3": (sib3, [("neighbour cells, same freq. (4)", "46, 69, 90, 100")], "20261126_test_si")},
            "plan": {"messages": [{"sibs": ["SIB3"], "T": 64, "frames": [0, 1]}], "window_slots": 20, "search_space": 0},
            "files": ["a.json"], "error": None}
    app.on_sibinfo({"key": (317, 2667750), "info": info})
    root.update()
    tops = [app.sib_tree.item(i, "text") for i in app.sib_tree.get_children()]
    assert tops == ["SIB1", "SIB3"], tops
    kids = [app.sib_tree.item(i, "text") for i in app.sib_tree.get_children(app.sib_tree.get_children()[0])]
    assert kids == ["network (PLMN) 1", "bands", "SI window"], kids
    assert app.nb_pcis == [46, 69, 90, 100] and app.nb_from == "SIB3 of PCI 317", (app.nb_pcis, app.nb_from)
    assert app.nb_job() == {"pcis": [46, 69, 90, 100], "window_us": 15.0, "cancel": True}
    app.update_mib_tab()
    assert "SIB1 decoded" in app.mib_notes.get(), app.mib_notes.get()
    # neighbour curves: PCI 100 peaks 20 samples after the centre (1 us at 20 MS/s), PCI 46 is noise
    rng = np.random.default_rng(1)
    for _ in range(300):
        c100 = 1 + 0.1 * rng.random(601)
        c100[320] += 10
        app.on_result(dict(pci=317, ok=False, decode_ms=5.0, run=app.run_id, cfo=13.9e3, snr_db=-3.0, t_pss_us=None,
                           nb={"serving": 317, "fs": 20e6, "curves": {100: (c100, None), 46: (1 + 0.1 * rng.random(601),
                                                                                            None)}}))
    app.nb_last_draw = 0
    app.update_nb_tab()
    rows = {app.nb_tree.item(i, "values")[0]: app.nb_tree.item(i, "values") for i in app.nb_tree.get_children()}
    assert rows["100"][5] == "detected" and rows["100"][2] == "+1.000", rows["100"]
    assert rows["46"][5] == "" and rows["46"][2] == "", rows["46"]
    assert "relative to PCI 317: 300 cuts" in app.nb_head.get(), app.nb_head.get()
    app.close()
    print("tabs with injected data ok")


def replay(shots=None):
    root = tk.Tk()
    app = App(root, make_args(replay=str(SESSION), ssb=None))
    app.start()
    # SIB info of PCI 317 from results/ fills the list (15 PCIs from its SIB3), then the sums grow
    assert pump(root, lambda: len(app.nb_pcis) == 15, 90), ("no SIB3 neighbours", app.nb_pcis, app.sib_head.get())
    assert pump(root, lambda: app.nb_acc.get(317) is not None and app.nb_acc[317].n >= 1500, 400), "too few cuts"
    rows = {r["pci"]: r for r in app.nb_acc[317].rows(app.nb_fs)}
    want = {100: 1.875, 90: 0.567, 69: 0.118}
    for p, tau in want.items():
        assert rows[p]["stat_sss"] >= 7, (p, rows[p]["stat_sss"])
        assert abs(rows[p]["tau_us"] - tau) < 0.03, (p, rows[p]["tau_us"])
    print("replay: " + ", ".join(f"PCI {p} {rows[p]['tau_us']:+.3f} us (SSS {rows[p]['stat_sss']:.1f})" for p in want)
          + f" from {app.nb_acc[317].n} cuts")
    if shots:
        app.nb_last_draw = 0
        app.update_nb_tab()
        app.tabs.select(app.tab_nb)
        grab(root, Path(shots) / "gui_neighbours.png")
        app.tabs.select(app.tab_sib)
        grab(root, Path(shots) / "gui_system_info.png")
    app.close()


def main():
    injected()
    if SESSION.exists():
        replay(sys.argv[1] if len(sys.argv) > 1 else None)
    else:
        print("private recording not present: replay part skipped")


if __name__ == "__main__":   # the decoder processes (spawn) import this file
    main()
