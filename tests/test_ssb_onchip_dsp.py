"""The on-chip SSB detector (firmware ssb_dsp.h), compiled for the host, on recorded bursts:
PSS position, SSS/PCI on the chip, the cell search over all N_ID2 and CFO hypotheses, and the
decimated 5 MS/s cuts decoded by the NumPy receiver.

Reference: the MATLAB detector's results for the same captures (results/*_ssb.csv, or the
excerpt in samples/ with its expected_ssb.csv):
every SSB of PCI 165 that MATLAB decoded must be found, at the same position (to 0.3 us).

Run from the project root:  .venv_local/bin/python -m unittest tests.test_ssb_onchip_dsp -v
Skipped when the captures or a C compiler are missing.
"""
import csv
import json
import os
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


def write_input(tmp: Path, folder: Path):
    recs = [json.loads(line) for line in (folder / "bursts.jsonl").read_text().splitlines() if line.strip()]
    blob = bytearray()
    with open(folder / "bursts.bin", "rb") as f:
        for r in recs:
            f.seek(r["offset"])
            w = raw_words(f.read(r["nbytes"]), r["n"], r["fmt"])
            blob += np.uint32(len(w)).tobytes() + w.astype("<u4").tobytes()
    inp = tmp / f"{folder.name}.bin"
    inp.write_bytes(bytes(blob))
    return recs, inp


def harness(exe: Path, inp: Path, fs: float, *args, exact: bool = False):
    env = dict(os.environ, SSB_SEARCH="exact" if exact else "fft")
    out = subprocess.run([str(exe), str(inp), str(int(fs))] + [str(a) for a in args], check=True,
                         capture_output=True, text=True, env=env).stdout
    # idx n_out lag stat cfo raw_idx ddc_ms corr_ms nid2 cfo_off pci sss_stat sss_ratio sss_ms
    return [line.split() for line in out.splitlines()]


def run_session(exe: Path, tmp: Path, folder: Path, fs: float, shift: int, nid2: int):
    recs, inp = write_input(tmp, folder)
    rows = {}
    for p in harness(exe, inp, fs, shift, nid2):
        rows[recs[int(p[0])]["idx"]] = dict(stat=float(p[3]), cfo=float(p[4]), raw=int(p[5]), ddc_ms=float(p[6]),
                                            corr_ms=float(p[7]), pci=int(p[10]), sss=float(p[11]))
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
                    self.assertEqual(got["pci"], 165, f"burst {b}: PCI from the on-chip SSS")
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


    def _sample20(self):
        folder, csv_path = find_session("20261123_112009_n1_2131")
        if folder is None:
            self.skipTest("no capture of the 20 MS/s session (data/ or samples/)")
        exp = {int(r["burst"]): r for r in matlab_ssbs(csv_path, "165")}
        recs, inp = write_input(self.tmp, folder)
        return recs, inp, exp

    def test_cell_search(self):
        """SSBDET with N_ID2 3 and a CFO span: from 0 kHz, all N_ID2, -15..+15 kHz (FFT search, 3.9 kHz
        bins), against the time-domain search over 5 kHz steps."""
        recs, inp, exp = self._sample20()
        rows = [p for p in harness(self.exe, inp, 20e6, int(2125.45e6 - 2131e6), 3, 15000) if float(p[3]) >= THRESHOLD]
        found = {recs[int(p[0])]["idx"]: p for p in rows}
        self.assertEqual(set(found), set(exp), "cell search: detections differ from MATLAB's SSBs")
        for b, p in found.items():
            self.assertEqual(int(p[8]), 0)
            self.assertEqual(int(p[10]), 165)
            # single-SSB PSS-halves estimate: noisy at low SNR; tracking needs it within about +-2.5 kHz
            self.assertLess(abs(int(p[9]) + float(p[4]) - 10840), 3000, "CFO of the cell search")
        cfos = [int(p[9]) + float(p[4]) for p in found.values()]
        self.assertLess(abs(np.median(cfos) - 10840), 800, "median CFO of the cell search")
        ref = {recs[int(p[0])]["idx"]: p for p in harness(self.exe, inp, 20e6, int(2125.45e6 - 2131e6), 3, 15000, exact=True)
               if float(p[3]) >= THRESHOLD}
        self.assertEqual(set(ref), set(found), "FFT and time-domain search find different SSBs")
        for b in found:
            self.assertLessEqual(abs(int(found[b][2]) - int(ref[b][2])), 1, f"burst {b}: lag of the FFT search")
        print(f"\ncell search: {len(found)}/{len(exp)} SSBs, N_ID2 0, PCI 165, "
              f"CFO {np.median([int(p[9]) + float(p[4]) for p in found.values()]):+.0f} Hz")

    def test_decimated_cuts(self):
        """The 5 MS/s cuts of the cut DDC decode in the NumPy receiver like the raw bursts."""
        from capture import nr_ssb as nr
        recs, inp, exp = self._sample20()
        cutfile = self.tmp / "cuts.bin"
        harness(self.exe, inp, 20e6, int(2125.45e6 - 2131e6 + 11200), 0, 0, cutfile)
        buf, pos, ok, tdiff, n_cut = cutfile.read_bytes(), 0, 0, [], 0
        while pos < len(buf):
            idx, n, start, R, pss, _ = np.frombuffer(buf[pos:pos + 24], dtype="<u4")
            pos += 24
            iq = np.frombuffer(buf[pos:pos + 4 * int(n)], dtype="<i2").astype(float)
            pos += 4 * int(n)
            start = int(np.int32(start))
            rec = recs[int(idx)]
            r = nr.receive(iq[0::2] + 1j * iq[1::2], 20e6 / int(R), (int(pss) - start) / int(R), 2125.45e6, 2125.45e6,
                           cfo=0.0, nid2=0, t_first_us=rec["t0_us"] + start / 20e6 * 1e6)
            n_cut += 1
            if r["ok"] and r["pci"] == 165:
                ok += 1
                if rec["idx"] in exp:
                    tdiff.append(r["t_pss_us"] - float(exp[rec["idx"]]["tPss_us"]))
        print(f"\ndecimated cuts: {n_cut} cuts at 5 MS/s, MIB of PCI 165 in {ok}; PSS time vs MATLAB "
              f"{np.median(tdiff):+.3f} us (spread {np.std(tdiff):.3f})")
        self.assertGreaterEqual(ok, len(exp) - 1)
        self.assertLess(abs(np.median(tdiff)), 0.1)


if __name__ == "__main__":
    unittest.main()
