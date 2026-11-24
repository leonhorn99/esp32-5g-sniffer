"""Channel selection in the live view, with the fake board (run on a virtual display by
tests/test_ssb_live.py): Track from the entries, a scan from the window, double-click on a
found cell, and that results of the previous channel are dropped."""
import time
import tkinter as tk
import types

import numpy as np

from capture.ssb_live import App, ScanSource, lo_for_ssb, nr_band, parse_los
from tests.fake_device import FakeDevice


def check_helpers():
    assert parse_los("2131 2670") == [2131, 2670]
    assert parse_los("2131, 2140-2150/5") == [2131, 2140, 2145, 2150]
    assert nr_band(2125.45) == "n1" and nr_band(2667.75) == "n41/n7"
    assert lo_for_ssb(2125.45, 2131) == 2131             # within 2..6.1 MHz: kept
    assert lo_for_ssb(2667.75, 2131) == 2673             # elsewhere: about 5 MHz above the SSB
    assert lo_for_ssb(2727.0, None) == 2722              # near the top of the range: below


def main():
    check_helpers()
    dev = FakeDevice(caps=("RXLIMITS", "TUNEEXT", "BIGCAP", "CAPT", "FWOPT", "SSBDET", "SSBDET3", "BOOT"),
                     junk_before_sync=False)
    dev.start()
    args = types.SimpleNamespace(
        replay=None, port=dev.port, freq=2131, ssb=2125.45, nid2=None, cfo=None, gain="40", rate=20, bits=8, thr=14.0,
        pre=20, post=20, period=20.0, every=1, scan_bursts=150, cfo_span=30000, raw_cuts=False, record=False, label="t",
        note="", workers=1, bandwidth=None, dccal=False, speed=1.0, hide=[], geometry="1560x920", dpi=90,
        autostart=False, quit_after=None, snapshot=None)
    root = tk.Tk()
    app = App(root, args)
    app.scan_bursts = 2


    def pump(cond, timeout=20.0):
        t_end = time.monotonic() + timeout
        while time.monotonic() < t_end:
            root.update()
            if cond():
                return True
            time.sleep(0.02)
        return False


    def ssbdet_cmds():
        return [line.split() for line in dev.log if line.startswith("SSBDET")]


    # 1. a channel from a scan: LO 2670, SSB 2667.75, N_ID2 2, CFO 13.93 kHz -> chip searches +-10 kHz
    app.track(2670, 2667.75, nid2=2, cfo=13930.0, what="PCI 317")
    assert pump(lambda: app.n_frames >= 5), "no frames after Track"
    assert app.cfg["f_lo"] == 2670e6 and app.cfg["f_ssb"] == 2667.75e6, app.cfg
    assert "FREQ 2670" in dev.log
    c = ssbdet_cmds()[-1]                 # SSBDET rate offset nid2 cfo thr pre post period bits every max track dec span
    assert (int(c[2]), int(c[3]), int(c[4]), int(c[14])) == (-2250000, 2, 13930, 10000), c
    assert int(c[13]) == 1, "decimated cuts with SSBDET3"
    run_before = app.run_id

    # 2. a scan from the window (fake bursts hold no SSB): it stops the stream, scans, ends cleanly
    app.v_scan.set("2131 3500")           # 3500 MHz is not tunable: skipped
    app.scan_entry()
    assert pump(lambda: isinstance(app.source, ScanSource) and not app.source.is_alive()), "scan did not finish"
    assert pump(lambda: str(app.btn_start.cget("state")) == "normal")
    assert "FREQ 2131" in dev.log and not any(line == "FREQ 3500" for line in dev.log)
    log = app.log.get("1.0", "end")
    assert "skipping [3500]" in log and "scan LO 2131 MHz: no cell" in log and "scan finished" in log, log

    # 3. a found cell in the table; double-click tracks it (N_ID2 from the PCI, CFO from the scan)
    cell = {"f_ssb_mhz": 2125.45, "gscn": 5314, "pci": 166, "ssbs": 5, "mibs": 4, "snr_db": -4.1, "rsrp_dbfs": -58.0,
            "cfo_hz": 10990.0, "ssb_index": [0], "mib": None, "sfn": 10}
    app.on_scan({"kind": "scan", "lo": 2131, "gain": 40, "cells": [cell]})
    root.update()
    iid = "2131_2125.450_166"
    assert app.found.exists(iid) and app.found.item(iid)["values"][0] == "n1"
    mib = dict(sfn=1, hrf=0, scs_common_khz=15, kssb=6, dmrs_typeA_pos=3, pdcch_config_sib1=16, coreset0=1,
               searchspace0=0, cell_barred=False, intra_freq_reselection="allowed")
    app.on_result(dict(pci=317, ok=True, cfo=13.9e3, snr_db=4.0, mib=mib, issb=3, t_pss_us=2e6, t_frame_us=2e6 - 1000,
                       grid=np.ones((240, 4)), nvar=1.0, eq=np.ones(432, complex), channel=np.ones(240), decode_ms=5.0,
                       stat=30.0, sss_stat=50.0, epoch=0, env_t=None, run=app.run_id))
    app.redraw()
    root.update()
    assert app.cells.exists("317") and 317 in app.cfo_pts
    app.on_found_double(types.SimpleNamespace(y=app.found.bbox(iid)[1] + 5))
    assert pump(lambda: app.n_frames >= 5 and app.cfg.get("f_lo") == 2131e6), "no frames after the double-click"
    c = ssbdet_cmds()[-1]
    assert (int(c[2]), int(c[3]), int(c[4])) == (-5550000, 1, 10990), c
    assert app.run_id > run_before
    assert not app.cells.exists("317") and 317 not in app.cfo_pts, "the previous channel's PCI is still shown"

    # 4. a decode result of the earlier channel is dropped
    n = app.n_decoded
    app.on_result({"run": run_before, "pci": 317, "decode_ms": 5.0})
    assert app.n_decoded == n and 317 not in app.series

    # 4b. right-click "track on the chip" for a PCI seen through the PC's extra N_ID2 tries
    app.series[165] = {"t": [1.0], "cfo": [12.1], "snr": [-1.0], "ok": [True], "tf": [], "sfn": [], "tfr": []}
    app.track_cell(165)
    assert pump(lambda: app.n_frames >= 5 and ssbdet_cmds()[-1][3] == "0"), "track_cell"
    c = ssbdet_cmds()[-1]
    assert (int(c[2]), int(c[3]), int(c[4]), int(c[14])) == (-5550000, 0, 12100, 10000), c

    # 5. manual entry: SSB outside the LO's span -> LO moved next to it
    app.v_lo.set("2131")
    app.v_ssb.set("2667.75")
    app.track_entry()
    assert pump(lambda: app.cfg.get("f_lo") == 2673e6), "manual Track"
    assert "using LO 2673 MHz" in app.log.get("1.0", "end") or app.args.freq == 2673

    app.close()
    dev.stop()
    print("channel selection ok")


if __name__ == "__main__":   # the decoder processes (spawn) import this file
    main()
