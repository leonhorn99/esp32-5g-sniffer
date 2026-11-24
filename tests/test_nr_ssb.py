"""The NumPy SSB receiver (capture/nr_ssb.py) against the 5G Toolbox and the MATLAB detector.

Run from the project root:  .venv_local/bin/python -m unittest tests.test_nr_ssb -v
Reference vectors: tests/fixtures/nr_pbch_vectors.json (matlab/export_pbch_tables.m).
The session test uses data/20261123_140216_ssbdet_v2 with its MATLAB result, or the excerpt in samples/.
"""
import csv
import json
import time
import unittest
from pathlib import Path

import numpy as np

from capture import nr_ssb as nr
from tests.data_paths import find_session

ROOT = Path(__file__).resolve().parents[1]
VEC = json.loads((ROOT / "tests" / "fixtures" / "nr_pbch_vectors.json").read_text())


def cplx(a):
    a = np.asarray(a, float)
    return a[0] + 1j * a[1] if a.shape[0] == 2 else a[:, 0] + 1j * a[:, 1]


class TestAgainstToolbox(unittest.TestCase):
    def test_sequences_and_indices(self):
        for v in VEC:
            pci = v["pci"]
            np.testing.assert_array_equal(nr.pss(pci % 3), v["pss"])
            np.testing.assert_array_equal(nr.sss(pci), v["sss"])
            ibar = v["issb"] + 4 * v["hrf"]
            np.testing.assert_allclose(nr.pbch_dmrs(pci, ibar), cplx(v["dmrs"]), atol=1e-12)
            pbch_i, dmrs_i = nr.ssb_indices(pci % 4)
            np.testing.assert_array_equal(pbch_i, v["pbchInd"])
            np.testing.assert_array_equal(dmrs_i, v["dmrsInd"])

    def test_bch_and_pbch_encoding(self):
        for v in VEC:
            cw = nr.bch_encode(v["trblk"], v["sfn"], v["hrf"], v["kssb"], v["pci"])
            np.testing.assert_array_equal(cw, v["cw"], err_msg=f"nrBCH, PCI {v['pci']}")
            sym = nr.qpsk(cw ^ nr.pbch_scrambling(v["pci"], v["issb"]))
            np.testing.assert_allclose(sym, cplx(v["pbch"]), atol=1e-12, err_msg="nrPBCH")

    def test_bch_decoding_clean(self):
        for v in VEC:
            llr = 1.0 - 2.0 * np.asarray(v["cw"], float)
            mib = nr.bch_decode(llr, v["pci"])
            self.assertIsNotNone(mib)
            self.assertEqual(mib["bits"], "".join(str(b) for b in v["trblk"]))
            self.assertEqual(mib["sfn"] & 15, int("".join(str(b) for b in v["dec"]["lsbsfn"]), 2))
            self.assertEqual(mib["sfn"], v["sfn"] & 15 | (int(mib["bits"][1:7], 2) << 4))
            self.assertEqual(mib["hrf"], v["dec"]["hrf"])
            self.assertEqual(mib["kssb"] >> 4, v["dec"]["msbidxoffset"] >> 4 if v["dec"]["msbidxoffset"] > 1
                             else v["dec"]["msbidxoffset"])

    def test_bch_decoding_noisy_and_speed(self):
        rng = np.random.default_rng(3)
        v = VEC[2]
        x = 1.0 - 2.0 * np.asarray(v["cw"], float)
        ok, t = 0, []
        for _ in range(20):
            sigma = 10 ** (2 / 20)           # Es/N0 per coded bit -2 dB (the MIB is rate 56/864)
            llr = 2 * (x + sigma * rng.standard_normal(864)) / sigma ** 2
            t0 = time.perf_counter()
            mib = nr.bch_decode(llr, v["pci"])
            t.append(time.perf_counter() - t0)
            ok += mib is not None and mib["bits"] == "".join(str(b) for b in v["trblk"])
        print(f"\nBCH list decoding at -2 dB per coded bit: {ok}/20, median {1e3 * np.median(t):.1f} ms")
        self.assertGreaterEqual(ok, 19)

    def test_crc_rejects_noise(self):
        # (a wrong PCI would not be caught here: the payload scrambling comes before the CRC;
        # in the receiver the PCI is fixed by the SSS and the PBCH scrambling before decoding)
        rng = np.random.default_rng(5)
        for _ in range(20):
            self.assertIsNone(nr.bch_decode(rng.standard_normal(864), 165))

    def test_symbol_times(self):
        self.assertAlmostEqual(nr.useful_start_s(2), 4544 / 30.72e6)
        self.assertAlmostEqual(nr.useful_start_s(16), 1e-3 + 4544 / 30.72e6)


class TestScs30(unittest.TestCase):
    """30 kHz SSBs, case C, L_max 8 (TDD bands): a burst from nrWaveformGenerator with the RF
    symbol phase (matlab/export_scs30_fixture.m)."""

    def test_case_c_lmax8(self):
        fx = json.loads((ROOT / "tests" / "fixtures" / "ssb_caseC_30k.json").read_text())
        v = np.fromfile(ROOT / "tests" / "fixtures" / "ssb_caseC_30k.bin", dtype="<i2").astype(float)
        x, fs = v[0::2] + 1j * v[1::2], fx["fs"]
        self.assertAlmostEqual(nr.useful_start_s(2, 30e3), fx["ssb"][0]["t_useful_s"])
        peaks = nr.pss_search(x, fs, 0.0, cfo_max=0, scs=30e3, max_peaks=10)
        self.assertEqual(len(peaks), 8)
        for p, s in zip(peaks, fx["ssb"]):
            self.assertLess(abs(p["t_pss"] / fs - s["t_useful_s"]), 1 / fs)
            r = nr.receive(x, fs, p["t_pss"], fx["f_ssb"], fx["f_ssb"], cfo=0, t_first_us=0.0, scs=30e3, lmax=8,
                           pattern="C")
            self.assertTrue(r["ok"])
            self.assertEqual((r["pci"], r["issb"], r["mib"]["sfn"], r["mib"]["hrf"]), (fx["pci"], s["issb"], s["sfn"], s["hrf"]))
            self.assertLess(abs(r["cfo"]), 20)
            self.assertLess(abs(r["t_frame_us"]), 0.05)


SESSION = "20261123_140216_ssbdet_v2"
SESSION_DIR, EXPECTED = find_session(SESSION)


@unittest.skipUnless(SESSION_DIR is not None, "session (or its sample) not present")
class TestOnRecordedCuts(unittest.TestCase):
    """The first 150 cuts of the SSBDET v2 session against MATLAB's ssb_detect (with the
    symbol-phase fix): same PCI and SFN, CFO and PSS time close."""

    def test_session(self):
        from capture.session import load_session
        meta, recs, iq = load_session(SESSION_DIR)
        ref = {}
        with open(EXPECTED) as fh:
            for r in csv.DictReader(fh):
                if r["sssOk"] == "1":
                    ref.setdefault(int(r["burst"]), []).append(r)
        f_lo = meta["settings"]["freq_mhz"][0] * 1e6
        f_ssb = meta["settings"]["ssb_mhz"] * 1e6
        n_cmp = n_ok_py = n_ok_m = both = 0
        cfo_d, t_d, times = [], [], []
        for rec, x in list(zip(recs, iq))[:150]:
            t0 = time.perf_counter()
            out = nr.receive(x, rec["fs"], rec["ssb_pss_in_cut"], f_ssb, f_lo, cfo=rec["ssb_cfo_hz"],
                             nid2=rec["ssb_nid2"], t_first_us=rec["t0_us"])
            times.append(time.perf_counter() - t0)
            m = [r for r in ref.get(rec["idx"], []) if r["nid2"] == str(rec["ssb_nid2"])]
            if not m:
                continue
            m = min(m, key=lambda r: abs(float(r["tPss_us"]) - out.get("t_pss_us", 0)))
            n_cmp += 1
            n_ok_py += out["ok"]
            n_ok_m += m["mibOk"] == "1"
            if out["ok"] and m["mibOk"] == "1":
                both += 1
                self.assertEqual(out["pci"], int(m["pci"]))
                self.assertEqual(out["mib"]["sfn"], int(m["sfn"]))
                self.assertEqual(out["issb"], int(m["iSSB"]))
                cfo_d.append(out["cfo"] - float(m["cfo_Hz"]))
                t_d.append(out["t_pss_us"] - float(m["tPss_us"]))
        print(f"\n{SESSION}, first 150 cuts: MIB decoded by Python {n_ok_py}, by MATLAB {n_ok_m} "
              f"(of {n_cmp} SSBs MATLAB found), both {both}; CFO Python - MATLAB median {np.median(cfo_d):+.0f} Hz "
              f"(spread {np.std(cfo_d):.0f}); PSS time Python - MATLAB median {np.median(t_d):+.3f} us "
              f"(spread {np.std(t_d):.3f}); {1e3 * np.median(times):.1f} ms per cut")
        self.assertGreaterEqual(n_ok_py, 0.97 * n_ok_m)
        self.assertLess(abs(np.median(cfo_d)), 150)
        self.assertLess(abs(np.median(t_d)), 0.1)


if __name__ == "__main__":
    unittest.main()
