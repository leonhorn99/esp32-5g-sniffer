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


class _Clock:
    """Virtual time for capture_lo: its sleeps advance this clock instead of waiting."""
    def __init__(self):
        self.t = 0.0

    def sleep(self, s):
        self.t += s


class _Board:
    """Bursts of 32704 samples at 20 MS/s; each capture takes `spacing_ms` (+- jitter) of the clock."""
    max_samples = 32704

    def __init__(self, clock, spacing_ms=94.0, jitter_ms=1.0, seed=1):
        self.clock, self.spacing, self.jitter = clock, spacing_ms, jitter_ms
        self.rng = np.random.default_rng(seed)
        self.t0s = []

    def capture(self, n, rate_code, bits=8):
        import types
        t0 = self.clock.t * 1e6
        self.t0s.append(t0)
        self.clock.t += (self.spacing + self.rng.uniform(-self.jitter, self.jitter)) / 1e3
        return types.SimpleNamespace(t0_us=t0, crc_ok=True, fs=20e6, fmt="CAP16", iq=lambda: np.zeros(n, complex))


class _Pool:
    """Runs a submitted job at once (the jobs here are stand-ins for scan_burst)."""
    def submit(self, fn, *a):
        from concurrent.futures import Future
        f = Future()
        f.set_result(fn(*a))
        return f


class TestSweep(unittest.TestCase):
    def test_coverage(self):
        from capture.cell_scan import Coverage, SSB_SPAN_US
        c = Coverage()
        c.add(0.0, 1635.2)                     # holds SSBs starting 0 .. 1305 us
        self.assertEqual(c.looks(), 0)
        self.assertEqual(int(c.seen.sum()), int((1635.2 - SSB_SPAN_US) // 10))
        for k in range(16):
            c.add(k * 1250.0, 1635.2)
        self.assertGreaterEqual(c.looks(), 1)
        c2 = Coverage()
        c2.add(0.0, 1635.2)
        s = c2.next_start(1, 1635.2)           # the gap is 1310 .. 20000 us: a burst centred on it
        self.assertAlmostEqual((s + (1635.2 - SSB_SPAN_US) / 2) % 20000, (1310 + 20000) / 2, delta=20)

    def test_timed_bursts(self):
        """Bursts timed into the gaps see every SSB position twice in about 60 bursts; with a spacing
        that is a multiple of the SSB period, untimed bursts would never get there."""
        from unittest import mock
        from capture import cell_scan as cs
        for spacing, jitter, most in ((94.0, 1.0, 72), (80.0, 0.0, 60)):
            clock = _Clock()
            board = _Board(clock, spacing, jitter)
            with mock.patch.object(cs.time, "sleep", clock.sleep), \
                    mock.patch.object(cs, "scan_burst", lambda *a: []):
                futs, n = cs.capture_lo(board, 2130, 2, 200, _Pool(), looks=2, say=lambda _t: None)
            cov = cs.Coverage()
            for t0 in board.t0s:
                cov.add(t0, 32704 / 20e6 * 1e6)
            self.assertGreaterEqual(cov.looks(), 2)
            self.assertLessEqual(n, most, f"spacing {spacing} ms")
            self.assertEqual(len(futs), n)
        clock = _Clock()
        board = _Board(clock, 80.0, 0.0)
        with mock.patch.object(cs.time, "sleep", clock.sleep), mock.patch.object(cs, "scan_burst", lambda *a: []):
            _, n = cs.capture_lo(board, 2130, 2, 100, _Pool(), say=lambda _t: None)      # untimed: all 100
        self.assertEqual(n, 100)

    def test_cases_and_cfo_window(self):
        from capture.cell_scan import cfo_window, ssb_cases
        from capture.ssb_sources import parse_los
        self.assertEqual(ssb_cases(2131), [(15e3, 4)])
        self.assertEqual(ssb_cases(2670), [(15e3, 4), (30e3, 8)])
        self.assertEqual(ssb_cases(5500), [(30e3, 8)])
        self.assertEqual(cfo_window(2131, None), (0.0, 30e3))
        self.assertEqual(cfo_window(5990, None), (0.0, 12 * 5990))
        c, s = cfo_window(2125.45, -5.8)
        self.assertAlmostEqual(c, 12327.6, places=1)
        self.assertEqual(s, 5e3)
        los = parse_los("2g 5g")
        self.assertEqual((len(los), los[0], los[50], los[51], los[-1]), (152, 2130, 2730, 4790, 5990))
        self.assertEqual(parse_los("2175-2199"), [2175, 2187, 2199])

    def test_sweep_learns_the_crystal(self):
        """A decoded cell at the first LO sets the crystal; the next LOs search +-2 ppm around it."""
        import types
        from unittest import mock
        from capture import cell_scan as cs
        seen = []

        def fake_scan(x, fs, f_lo, c0, span, t0, scale, nid2s, scs, lmax):
            seen.append((f_lo / 1e6, c0, span, scs))
            if f_lo == 2130e6:                 # PCI 165 at 2125.45 MHz, CFO +12.33 kHz: crystal -5.8 ppm
                return [{"f_ssb": 2125.45e6, "scs": scs, "pci": 165, "ok": True, "snr_db": 3.0, "rsrp_dbfs": -60.0,
                         "cfo": 12327.6, "issb": 0, "pss_stat": 40.0, "sss_stat": 60.0, "t_pss_us": 1.0, "sfn": 5,
                         "mib": {}}]
            return []
        clock = _Clock()
        args = types.SimpleNamespace(gain="40", bits=8)
        with mock.patch.object(cs.time, "sleep", clock.sleep), mock.patch.object(cs, "scan_burst", fake_scan), \
                mock.patch("capture.capture_bursts._setup", lambda *a: None):
            out = list(cs.sweep(_Board(clock), [2130, 2670, 4790], 2, _Pool(), args, looks=2, say=lambda _t: None))
        self.assertEqual([r["lo"] for r in out], [2130, 2670, 4790])
        self.assertEqual(out[0]["cells"][0]["pci"], 165)
        self.assertAlmostEqual(out[0]["ppm"], -5.8, places=3)
        self.assertEqual(out[0]["cfo"], (0.0, 30e3))                         # before: unknown, +-30 kHz
        c, s = out[1]["cfo"]
        self.assertAlmostEqual(c, 5.8 * 2670, places=0)                     # after: -ppm x f +- 2 ppm
        self.assertAlmostEqual(s, 2 * 2670)
        self.assertEqual({x[3] for x in seen if x[0] == 2670}, {15e3, 30e3})
        self.assertEqual({x[3] for x in seen if x[0] == 4790}, {30e3})
        self.assertTrue(all(r["bursts"] <= 80 for r in out))


class TestMibWords(unittest.TestCase):
    """capture/nr_mib.py on the MIBs of 24 Nov (TS 38.213 Tables 13-1 and 13-11)."""
    M165 = dict(sfn=560, hrf=0, scs_common_khz=15, kssb=6, dmrs_typeA_pos=3, pdcch_config_sib1=117, coreset0=7,
                searchspace0=5, cell_barred=False, intra_freq_reselection="allowed", bits="0" * 24, path_metric=1.0)

    def test_cells_of_24_nov(self):
        from capture.nr_mib import describe
        d = describe(self.M165, 2125.45e6, 0)
        self.assertEqual((d["coreset0"], d["coreset0_mhz"], d["ss0_slots"], d["ss0_sym"], d["sib1"]),
                         ("48 RB × 1 sym", "2120.68-2129.32", "5, 6 (even SFN)", "0", "yes"))
        d = describe(dict(self.M165, pdcch_config_sib1=16), 2667.75e6, 3)        # PCI 317, SSB index 3
        self.assertEqual((d["coreset0"], d["coreset0_mhz"], d["ss0_slots"]),
                         ("24 RB × 2 sym", "2665.50-2669.82", "3, 4 (even SFN)"))

    def test_occasions_and_edge_cases(self):
        from capture.nr_mib import content_key, describe, sib1_occasions
        self.assertEqual(sib1_occasions(self.M165, 1)["first_symbol"], 1)             # odd i: N_symb^CORESET
        self.assertEqual(sib1_occasions(dict(self.M165, pdcch_config_sib1=7 * 16 + 9), 3)["slots"], (1, 2))  # O 5, M 2
        self.assertEqual(sib1_occasions(dict(self.M165, pdcch_config_sib1=7 * 16 + 9), 3)["frame_parity"], 1)
        self.assertEqual(describe(dict(self.M165, kssb=30), 2125.45e6, 0)["coreset0"], "none")
        self.assertEqual(describe(dict(self.M165, pdcch_config_sib1=15 * 16), 2125.45e6, 0)["coreset0"], "reserved")
        self.assertEqual(describe(dict(self.M165, scs_common_khz=30), 2125.45e6, 0)["coreset0"], "not in table")
        self.assertEqual(content_key(self.M165), content_key(dict(self.M165, sfn=561, hrf=1, bits="1" * 24)))
        self.assertNotEqual(content_key(self.M165), content_key(dict(self.M165, cell_barred=True)))


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
