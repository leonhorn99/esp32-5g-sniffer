"""The on-chip SSB detector (firmware ssb_dsp.h), compiled for the host, on recorded bursts.

Reference: the MATLAB detector's results for the same captures (results/*_ssb.csv, or the
excerpt in samples/ with its expected_ssb.csv):
every SSB of PCI 165 that MATLAB decoded must be found, at the same position (to 0.3 us).

Run from the project root:  .venv_local/bin/python -m unittest tests.test_ssb_onchip_dsp -v
Skipped when the captures or a C compiler are missing.
"""
import csv
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

from capture.esp_sdr import raw_words
from tests.data_paths import find_session

ROOT = Path(__file__).resolve().parents[1]
FW = ROOT / "firmware" / "espsoup-firmware"
HARNESS = FW / "tests" / "ssb" / "host_ssb.c"
THRESHOLD = 14.0   # firmware default (thr_x10 = 140)

# session, fs, LO, SSB centre, CFO, N_ID2. The CFO is the within-symbol (PSS-halves) estimate;
# MATLAB's reported CFO (14854 / 16257 Hz) is biased by the gNB's per-symbol phase (TS 38.211 5.4).
CASES = [
    ("20261123_112009_n1_2131", 20e6, 2131e6, 2125.45e6, 11200, 0),
    ("20261123_113413_n1_2131_40M", 40e6, 2131e6, 2125.45e6, 12800, 0),
]


def build_harness(tmp: Path) -> Path:
    exe = tmp / "host_ssb"
    subprocess.run(["cc", "-O2", "-std=c11", "-Wall", "-Wextra", "-I", str(FW / "main/families/c5_c6_c61"),
                    str(HARNESS), "-lm", "-o", str(exe)], check=True, capture_output=True, text=True)
    return exe


def run_session(exe: Path, tmp: Path, folder: Path, fs: float, shift: int, nid2: int):
    session = folder.name
    recs = [json.loads(line) for line in (folder / "bursts.jsonl").read_text().splitlines() if line.strip()]
    blob = bytearray()
    with open(folder / "bursts.bin", "rb") as f:
        for r in recs:
            f.seek(r["offset"])
            w = raw_words(f.read(r["nbytes"]), r["n"], r["fmt"])
            blob += np.uint32(len(w)).tobytes() + w.astype("<u4").tobytes()
    inp = tmp / f"{session}.bin"
    inp.write_bytes(bytes(blob))
    out = subprocess.run([str(exe), str(inp), str(int(fs)), str(shift), str(nid2)], check=True,
                         capture_output=True, text=True).stdout
    rows = {}
    for line in out.splitlines():
        idx, n_out, lag, stat, cfo, raw_idx, ddc_ms, corr_ms = line.split()
        rows[recs[int(idx)]["idx"]] = dict(stat=float(stat), cfo=float(cfo), raw=int(raw_idx),
                                           ddc_ms=float(ddc_ms), corr_ms=float(corr_ms))
    t0 = {r["idx"]: r["t0_us"] for r in recs}
    return rows, t0


def matlab_ssbs(csv_path: Path, pci: str):
    with open(csv_path) as fh:
        return [r for r in csv.DictReader(fh) if r["sssOk"] == "1" and r["pci"] == pci]


@unittest.skipUnless(shutil.which("cc"), "no C compiler")
class TestOnChipDetector(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp())
        cls.exe = build_harness(cls.tmp)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_sessions(self):
        ran = 0
        for session, fs, f_lo, f_ssb, cfo, nid2 in CASES:
            folder, csv_path = find_session(session)
            if folder is None:
                continue
            ran += 1
            with self.subTest(session=session):
                shift = int(round(f_ssb - f_lo + cfo))
                rows, t0 = run_session(self.exe, self.tmp, folder, fs, shift, nid2)
                expected = matlab_ssbs(csv_path, "165")
                found, worst = 0, 0.0
                for r in expected:
                    b = int(r["burst"])
                    want = (float(r["tPss_us"]) - t0[b]) * fs / 1e6
                    got = rows[b]
                    err_us = (got["raw"] - want) / fs * 1e6
                    worst = max(worst, abs(err_us))
                    self.assertGreaterEqual(got["stat"], THRESHOLD, f"burst {b}: statistic {got['stat']:.1f}")
                    self.assertLess(abs(err_us), 0.3, f"burst {b}: PSS position off by {err_us:.2f} us")
                    found += 1
                exp_bursts = {int(r["burst"]) for r in expected}
                extra = sorted(b for b, v in rows.items() if v["stat"] >= THRESHOLD and b not in exp_bursts)
                stats = sorted(v["stat"] for b, v in rows.items() if b not in exp_bursts)
                ddc = np.median([v["ddc_ms"] for v in rows.values()])
                corr = np.median([v["corr_ms"] for v in rows.values()])
                print(f"\n{session}: {found}/{len(expected)} MATLAB SSBs found, worst position error {worst:.3f} us; "
                      f"detections in other bursts: {extra}; max statistic elsewhere {stats[-1]:.1f}; "
                      f"host time per burst: DDC {ddc:.2f} ms, PSS {corr:.2f} ms")
                self.assertLessEqual(len(extra), 3, f"too many detections outside MATLAB's SSBs: {extra}")
        if not ran:
            self.skipTest("no capture (data/ or samples/) present")


if __name__ == "__main__":
    unittest.main()
