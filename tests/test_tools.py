"""The raster PSS search, the cell scanner, the headless monitor and the SigMF export, on the
published samples (or the full captures when present).

Run from the project root:  .venv_local/bin/python -m unittest tests.test_tools -v
"""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from capture import nr_ssb as nr
from capture.session import load_session
from tests.data_paths import find_session

ROOT = Path(__file__).resolve().parents[1]
BURSTS, BURSTS_CSV = find_session("20261123_112009_n1_2131")
CUTS, _ = find_session("20261123_140216_ssbdet_v2")


def expected_bursts(csv_path):
    with open(csv_path) as fh:
        return {int(r["burst"]) for r in csv.DictReader(fh) if r["sssOk"] == "1" and r["pci"] == "165"}


class TestRaster(unittest.TestCase):
    def test_gscn(self):
        f = nr.gscn_raster(2125e6, 2126e6)
        self.assertIn(2125.45e6, f)
        self.assertTrue(all(abs((x / 50e3) - round(x / 50e3)) < 1e-9 for x in f))
        self.assertEqual(nr.gscn_raster(3000e6, 3002e6), [3000e6, 3001.44e6])
        from capture.cell_scan import gscn_number
        self.assertEqual(gscn_number(2125.45e6), 5314)
        self.assertEqual(gscn_number(2667.75e6), 6669)

    @unittest.skipUnless(BURSTS is not None, "capture not present")
    def test_raster_finds_every_ssb(self):
        meta, recs, iq = load_session(BURSTS)
        exp = expected_bursts(BURSTS_CSV)
        offs = [f - 2131e6 for f in nr.gscn_raster(2131e6 - 6.1e6, 2131e6 + 6.1e6)]
        found, false = set(), []
        for r, x in zip(recs, iq):
            for p in nr.pss_search_raster(x, 20e6, offs, np.arange(-5000, 25001, 5000)):
                o = nr.receive(x, 20e6, p["t_pss"], 2131e6 + p["f_off"], 2131e6, cfo=p["cfo"], nid2=p["nid2"])
                if o.get("pci") == 165 and o["ok"] and abs(p["f_off"] - (2125.45e6 - 2131e6)) < 1:
                    found.add(r["idx"])
                elif o.get("pci") is not None:
                    false.append((r["idx"], p["f_off"], o["pci"]))
        self.assertEqual(found & exp, exp & {r["idx"] for r in recs})
        self.assertEqual(false, [])


@unittest.skipUnless(BURSTS is not None, "capture not present")
class TestCellScan(unittest.TestCase):
    def test_session_scan(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "cells.json"
            subprocess.run([sys.executable, "-m", "capture.cell_scan", "--session", str(BURSTS), "--cfo", "10000",
                            "--out", str(out)], cwd=ROOT, check=True, capture_output=True, text=True, timeout=600)
            cells = json.loads(out.read_text())["los"][0]["cells"]
        c = cells[0]
        self.assertEqual((c["pci"], c["gscn"]), (165, 5314))
        self.assertAlmostEqual(c["f_ssb_mhz"], 2125.45)
        self.assertEqual(c["mibs"], len([b for b in expected_bursts(BURSTS_CSV)]))
        self.assertLess(abs(c["cfo_hz"] - 10840), 300)
        self.assertEqual(c["mib"]["pdcch_config_sib1"], 117)


@unittest.skipUnless(CUTS is not None, "capture not present")
class TestMonitorAndSigmf(unittest.TestCase):
    def test_monitor_replay(self):
        with tempfile.TemporaryDirectory() as d:
            jl = Path(d) / "m.jsonl"
            out = subprocess.run([sys.executable, "-m", "capture.ssb_monitor", "--replay", str(CUTS), "--speed", "0",
                                  "--jsonl", str(jl), "--expect-pci", "165", "--workers", "1"],
                                 cwd=ROOT, capture_output=True, text=True, timeout=600)
            rows = [json.loads(line) for line in jl.read_text().splitlines()]
        ssb = [r for r in rows if r["type"] == "ssb" and r["pci"] == 165 and r["ok"]]
        alerts = [r for r in rows if r["type"] == "alert"]
        self.assertGreaterEqual(len(ssb), 140, out.stderr[-500:])
        self.assertEqual({a["kind"] for a in alerts}, {"new_pci"}, alerts)        # only PCI 166, not expected
        self.assertTrue(all(a.get("pci") == 166 for a in alerts))
        resid = [r["timing_resid_us"] for r in ssb if r.get("timing_resid_us") is not None]
        self.assertLess(np.std(resid), 0.5)

    def test_sigmf(self):
        with tempfile.TemporaryDirectory() as d:
            base = Path(d) / "x"
            subprocess.run([sys.executable, "-m", "capture.sigmf_export", str(CUTS), str(base)], cwd=ROOT, check=True,
                           capture_output=True, text=True, timeout=300)
            meta = json.loads(base.with_suffix(".sigmf-meta").read_text())
            data = np.fromfile(base.with_suffix(".sigmf-data"), dtype=np.int8)
            _, recs, iq = load_session(CUTS)
        self.assertEqual(meta["global"]["core:datatype"], "ci8")
        self.assertEqual(len(meta["captures"]), len(recs))
        self.assertEqual(len(data) // 2, sum(len(x) for x in iq))
        first = data[:2 * len(iq[0])].astype(float)
        np.testing.assert_array_equal(first[0::2] + 1j * first[1::2], iq[0])
        self.assertEqual(len(meta["annotations"]), len(recs))


if __name__ == "__main__":
    unittest.main()
