"""capture/ssb_live.py: the decoder job, the clock fit, and a replay of the GUI on a virtual display.

Run from the project root:  .venv_local/bin/python -m unittest tests.test_ssb_live -v
The replay test needs xvfb-run; the session is data/20261123_140216_ssbdet_v2 or its excerpt in samples/.
"""
import os
import shutil
import subprocess
import sys
import unittest
from collections import deque
from pathlib import Path

import numpy as np

from capture.ssb_live import App, decode_job
from tests.data_paths import find_session

ROOT = Path(__file__).resolve().parents[1]
SESSION, _ = find_session("20261123_140216_ssbdet_v2")


class TestClockFit(unittest.TestCase):
    def test_quadratic_drift_removed(self):
        rng = np.random.default_rng(1)
        t = np.arange(0, 60e6, 20e3)                    # 60 s of SSBs every 20 ms (us, board clock)
        sfn = (np.arange(len(t)) * 2) % 1024
        ppm0, dppm = -7.2, -0.2 / 60e6                  # -7.2 ppm, changing by -0.2 ppm per minute
        clock = (ppm0 + dppm * t / 2) * 1e-6 * t      # integrated clock error
        tf = 1000.0 + t + clock + rng.uniform(-0.5, 0.5, len(t))   # 1 us timestamp rounding
        s = {"tf": deque(tf), "sfn": deque(sfn), "tfr": deque(t / 1e6)}
        # frame start on the board clock = true start (SFN x 10 ms) + clock error
        s["tf"] = deque(1000.0 + sfn * 1e4 + np.floor(np.arange(len(t)) * 2 / 1024) * 10.24e6 + clock
                        + rng.uniform(-0.5, 0.5, len(t)))
        _, r, ppm_now, ppm_min, rms = App.timing_residual(s)
        self.assertAlmostEqual(rms, 1 / np.sqrt(12), delta=0.02)
        self.assertAlmostEqual(ppm_now, ppm0 + dppm * t[-1], delta=0.01)
        self.assertAlmostEqual(ppm_min, -0.2, delta=0.01)


@unittest.skipUnless(SESSION is not None, "session not present")
class TestDecodeJob(unittest.TestCase):
    def test_cuts(self):
        from capture.session import load_session
        meta, recs, iq = load_session(SESSION)
        n_mib, pcis = 0, set()
        for rec, x in list(zip(recs, iq))[:20]:
            job = dict(iq=x, fs=rec["fs"], t_pss=rec["ssb_pss_in_cut"], nid2=rec["ssb_nid2"], cfo=rec["ssb_cfo_hz"],
                       f_ssb=2125.45e6, f_lo=2131e6, t_first_us=rec["t0_us"], seq=rec["seq"], mode="T", stat=20.0)
            for r in decode_job(job):
                pcis.add(r["pci"])
                n_mib += bool(r["ok"])
                self.assertEqual(r["grid"].shape, (240, 4))
                self.assertEqual(len(r["env_t"]), len(r["env_db"]))
        self.assertIn(165, pcis)
        self.assertGreaterEqual(n_mib, 18)


@unittest.skipUnless(SESSION is not None and shutil.which("xvfb-run"), "session or xvfb-run not present")
class TestGuiReplay(unittest.TestCase):
    def test_replay(self):
        out = subprocess.run(["xvfb-run", "-a", "-s", "-screen 0 1600x1000x24", sys.executable, "-m", "capture.ssb_live",
                              "--replay", str(SESSION), "--speed", "0", "--quit-after", "20"],
                             cwd=ROOT, capture_output=True, text=True, timeout=120)
        line = [x for x in out.stdout.splitlines() if x.startswith("ssb_live:")]
        self.assertTrue(line, out.stdout + out.stderr)
        print("\n" + line[0])
        self.assertIn("PCI 165", line[0])


@unittest.skipUnless(shutil.which("xvfb-run"), "xvfb-run not present")
class TestPciSelection(unittest.TestCase):
    def test_show_hide(self):
        out = subprocess.run(["xvfb-run", "-a", sys.executable, str(ROOT / "tests" / "gui" / "check_pci_selection.py")],
                             cwd=ROOT, capture_output=True, text=True, timeout=120,
                             env={**os.environ, "PYTHONPATH": str(ROOT)})
        self.assertIn("PCI selection OK", out.stdout, out.stdout + out.stderr)


@unittest.skipUnless(shutil.which("xvfb-run"), "xvfb-run not present")
class TestChannelSelection(unittest.TestCase):
    def test_track_scan_switch(self):
        """Track, a scan and a double-click on a found cell, with the fake board."""
        out = subprocess.run(["xvfb-run", "-a", sys.executable, str(ROOT / "tests" / "gui" / "check_channel_select.py")],
                             cwd=ROOT, capture_output=True, text=True, timeout=240,
                             env={**os.environ, "PYTHONPATH": str(ROOT)})
        self.assertIn("channel selection ok", out.stdout, out.stdout + out.stderr)


if __name__ == "__main__":
    unittest.main()
