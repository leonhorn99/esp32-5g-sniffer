"""Tests of the capture client against the protocol description and a fake board.

Run from the project root:  .venv_local/bin/python -m unittest discover -s tests -t . -v
"""
import struct
import tempfile
import unittest
import zlib
from pathlib import Path

import numpy as np

from capture.esp_sdr import EspSdr, pack_iq, payload_bytes, rate_code_for, unpack_iq
from capture.session import SessionWriter, load_session
from tests.fake_device import FakeDevice, pattern


class TestUnpack(unittest.TestCase):
    def test_cap16_layout(self):
        # Spec: byte 0 = I (int8), byte 1 = -Q (int8)
        payload = bytes([0x05, 0xFD, 0x80, 0x7F])  # (5, -3), (-128, 127)
        iq = unpack_iq(payload, 2, "CAP16")
        np.testing.assert_array_equal(iq.real, [5, -128])
        np.testing.assert_array_equal(iq.imag, [3, -127])  # Q = -(-Q)

    def test_cap20_layout(self):
        # Spec: 20-bit samples LSB first, I = bits 0-9, -Q = bits 10-19, signed 10 bit.
        s0 = (1 & 0x3FF) | ((-1 & 0x3FF) << 10)       # I=1,    -Q=-1
        s1 = (-512 & 0x3FF) | ((511 & 0x3FF) << 10)   # I=-512, -Q=511
        word = s0 | (s1 << 20)
        payload = word.to_bytes(5, "little")
        iq = unpack_iq(payload, 2, "CAP20")
        np.testing.assert_array_equal(iq.real, [1, -512])
        np.testing.assert_array_equal(iq.imag, [1, -511])

    def test_cap20_odd_length_bytes(self):
        self.assertEqual(payload_bytes(3, "CAP20"), 8)   # ceil(60/8)
        self.assertEqual(payload_bytes(32704, "CAP20"), 81760)
        self.assertEqual(payload_bytes(32704, "CAP16"), 65408)

    def test_roundtrip_random(self):
        rng = np.random.default_rng(1)
        for fmt, lim in (("CAP16", 128), ("CAP20", 512)):
            for n in (1, 2, 7, 1000, 1001):
                i = rng.integers(-lim, lim, n)
                q = rng.integers(-lim + 1, lim, n)  # -Q must fit, so Q > -lim
                iq = unpack_iq(pack_iq(i, q, fmt), n, fmt)
                np.testing.assert_array_equal(iq.real, i)
                np.testing.assert_array_equal(iq.imag, q)

    def test_rate_codes(self):
        self.assertEqual(rate_code_for(20e6), 2)
        self.assertEqual(rate_code_for(4e6), 5)


class TestFakeBoard(unittest.TestCase):
    def setUp(self):
        self.dev = FakeDevice(corrupt_seqs={2})
        self.dev.start()
        self.sdr = EspSdr(self.dev.port, timeout=2.0)

    def tearDown(self):
        self.sdr.ser.close()
        self.dev.stop()

    def test_connect_and_capture(self):
        info = self.sdr.connect()
        self.assertEqual(self.sdr.family, "C5SDR")
        self.assertEqual(self.sdr.max_samples, 32704)
        self.assertIn("CAPT", info["caps"])
        self.sdr.set_freq_mhz(2138)
        self.sdr.set_gain(50)
        b = self.sdr.capture(32704, 2, bits=10)   # seq 0
        self.assertTrue(b.crc_ok)
        self.assertEqual(b.t0_us, 1_000_000)
        i, q = pattern(32704, "CAP20", 0)
        iq = b.iq()
        np.testing.assert_array_equal(iq.real, i)
        np.testing.assert_array_equal(iq.imag, q)
        b = self.sdr.capture(1001, 5, bits=8)     # seq 1
        i, q = pattern(1001, "CAP16", 1)
        np.testing.assert_array_equal(b.iq().real, i)
        np.testing.assert_array_equal(b.iq().imag, q)
        b = self.sdr.capture(500, 2, bits=10)     # seq 2 is corrupted on purpose
        self.assertFalse(b.crc_ok)
        self.assertAlmostEqual(self.sdr.temperature(), 43.9)

    def test_without_capt_uses_cap(self):
        self.dev.caps = [c for c in self.dev.caps if c != "CAPT"]
        self.sdr.connect()
        b = self.sdr.capture(100, 2, bits=10)
        self.assertTrue(b.reply.startswith("DATA "))
        self.assertIsNone(b.t0_us)
        self.assertTrue(b.crc_ok)

    def test_session_roundtrip(self):
        device = self.sdr.connect()
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d) / "s"
            w = SessionWriter(folder, device, {"x": 1})
            for n, bits in ((101, 10), (64, 8), (500, 10)):
                w.add(self.sdr.capture(n, 2, bits=bits), freq_mhz=2138, gain=50, bandwidth_mhz=0)
            w.close()
            meta, recs, iq = load_session(folder)              # drops the corrupted burst
            self.assertEqual(len(recs), 2)
            self.assertEqual(meta["device"]["info"], "C5SDR 6 burst 16380")
            i, q = pattern(101, "CAP20", 0)
            np.testing.assert_array_equal(iq[0].real, i)
            np.testing.assert_array_equal(iq[0].imag, q)
            _, recs_all, _ = load_session(folder, only_crc_ok=False)
            self.assertEqual(len(recs_all), 3)


class TestCli(unittest.TestCase):
    def test_capture_cli(self):
        from capture.capture_bursts import main
        dev = FakeDevice()
        dev.start()
        try:
            with tempfile.TemporaryDirectory() as d:
                main(["capture", "--port", dev.port, "--freq", "2138", "2158", "--per-freq", "2",
                      "--rate", "20", "--bits", "10", "--gain", "50", "--count", "5",
                      "--out", d, "--label", "t", "--dccal"])
                folder = next(Path(d).iterdir())
                meta, recs, iq = load_session(folder)
                self.assertEqual(len(recs), 5)
                self.assertEqual([r["freq_mhz"] for r in recs], [2138, 2138, 2158, 2158, 2138])
                self.assertEqual(meta["settings"]["n"], 32704)
                self.assertIn("FREQ 2158", dev.log)
                self.assertIn("DCCAL", dev.log)
        finally:
            dev.stop()


    def test_capture_auto_gain(self):
        from capture.capture_bursts import main
        dev = FakeDevice(gain_scaling=True)   # signal = pattern * gain/80; pattern peak 512 counts
        dev.start()
        try:
            with tempfile.TemporaryDirectory() as d:
                main(["capture", "--port", dev.port, "--freq", "2142", "--rate", "20", "--bits", "10",
                      "--gain", "auto", "--count", "2", "--out", d])
                meta, recs, iq = load_session(next(Path(d).iterdir()))
                # highest index in steps of 4 with peak <= 255.5: 512*g/80 <= 255.5 -> g = 36
                self.assertEqual(meta["settings"]["gain"], 36)
                self.assertEqual(len(meta["settings"]["gain_scan"]), 21)
                self.assertTrue(all(r["gain"] == 36 for r in recs))
                self.assertLessEqual(np.max(np.abs(iq[0].real)), 256)
        finally:
            dev.stop()


class TestLineReader(unittest.TestCase):
    def test_partial_line_survives_timeout(self):
        """A header split by a pause must not be lost (cause of the 110 lost frames on 23 Nov)."""
        import os
        dev = FakeDevice()
        try:
            sdr = EspSdr(dev.port, timeout=2.0)
            os.write(dev.master, b"SSBF 1 0 000")
            with self.assertRaises(TimeoutError):
                sdr._read_line(timeout=0.2)
            os.write(dev.master, b"00000 5 0 0 20000000 16 T 150 0 11000 100\n")
            self.assertEqual(sdr._read_line(timeout=1.0), "SSBF 1 0 00000000 5 0 0 20000000 16 T 150 0 11000 100")
            sdr.ser.close()
        finally:
            dev.stop()


class TestSsbStream(unittest.TestCase):
    CAPS = ("RXLIMITS", "TUNEEXT", "GPIO", "BIGCAP", "CAPT", "TEMP", "FWOPT", "DCCAL", "SSBDET")

    def test_library(self):
        from capture.esp_sdr import ssb_stream
        dev = FakeDevice(caps=self.CAPS)
        dev.start()
        try:
            sdr = EspSdr(dev.port, timeout=2.0)
            sdr.connect()
            out = list(ssb_stream(sdr, 2, -5550000, 0, 11200, bits=8, max_frames=4))
            kinds = [k for k, _ in out]
            self.assertEqual(kinds, ["start", "frame", "frame", "frame", "frame", "end"])
            frames = [f for k, f in out if k == "frame"]
            self.assertEqual([f.n for f in frames], [600, 600, 0, 600])     # third is header-only
            self.assertTrue(all(f.crc_ok for f in frames))
            self.assertEqual(frames[0].mode, "A")
            self.assertAlmostEqual(frames[0].t_first_us, frames[0].t0_us + 1000 / 20e6 * 1e6)
            i, q = pattern(600, "CAP16", 1)
            np.testing.assert_array_equal(frames[1].iq().real, i)
            np.testing.assert_array_equal(frames[1].iq().imag, q)
            self.assertIn("SSBDET 2 -5550000 0 11200 140 20 20 20000 16 1 4 1", dev.log)
            sdr.ser.close()
        finally:
            dev.stop()

    def test_ssb3_fields_and_decimated_cuts(self):
        """Firmware ssb3: sub-us t0, on-chip PCI, decimated cuts, BOOT?; the recorder stores the cut rate."""
        from capture.esp_sdr import ssb_stream
        from capture.session import SsbRecorder, load_session
        dev = FakeDevice(caps=self.CAPS + ("SSBDET3", "BOOT"))
        dev.start()
        try:
            sdr = EspSdr(dev.port, timeout=2.0)
            device = sdr.connect()
            self.assertEqual(device["boot"], "BOOT 1 poweron 12345678")
            frames = [f for k, f in ssb_stream(sdr, 2, -5550000, 3, 11200, max_frames=4, dec=True, cfo_span=15000)
                      if k == "frame"]
            self.assertIn("SSBDET 2 -5550000 3 11200 140 20 20 20000 16 1 4 1 1 15000", dev.log)
            f = frames[0]
            self.assertEqual((f.pci, f.R, f.gsh, f.cfo_ddc_hz, f.t0_sub_ns), (165, 4, 2, 11000, 437))
            self.assertEqual(f.fs_cut, 5e6)
            self.assertAlmostEqual(f.t_first_us, f.t0_us + 0.437 + 1000 / 20e6 * 1e6)
            self.assertAlmostEqual(f.pss_in_cut, (1094 - 1000) / 4)
            with tempfile.TemporaryDirectory() as d:
                rec = SsbRecorder(Path(d) / "s", device, {}, 2131, 50, None, 2)
                for fr in frames:
                    rec.add(fr)
                rec.close()
                meta, recs, iq = load_session(Path(d) / "s")
                self.assertEqual(recs[0]["fs"], 5e6)
                self.assertEqual(recs[0]["ssb_pss_in_cut"], 23.5)
                self.assertEqual(recs[0]["ssb_pci_chip"], 165)
            sdr.ser.close()
        finally:
            dev.stop()

    def test_old_firmware_refuses_ssb3_options(self):
        from capture.esp_sdr import EspSdrError, ssb_stream
        dev = FakeDevice(caps=self.CAPS)
        dev.start()
        try:
            sdr = EspSdr(dev.port, timeout=2.0)
            sdr.connect()
            with self.assertRaises(EspSdrError):
                list(ssb_stream(sdr, 2, -5550000, 0, 11200, max_frames=2, dec=True))
            sdr.ser.close()
        finally:
            dev.stop()

    def test_cli_session(self):
        from capture.capture_bursts import main
        dev = FakeDevice(caps=self.CAPS)
        dev.start()
        try:
            with tempfile.TemporaryDirectory() as d:
                main(["ssb-stream", "--port", dev.port, "--freq", "2131", "--ssb", "2125.45", "--nid2", "0",
                      "--cfo", "11200", "--gain", "50", "--count", "6", "--duration", "10", "--out", d, "--label", "t"])
                folder = next(Path(d).iterdir())
                meta, recs, iq = load_session(folder)
                self.assertEqual(len(recs), 4)                     # 6 frames, 2 header-only
                self.assertEqual(meta["settings"]["offset_hz"], -5550000)
                self.assertAlmostEqual(recs[0]["t0_us"], recs[0]["capture_t0_us"] + 1000 / 20e6 * 1e6, places=3)
                self.assertEqual(recs[0]["ssb_pss_in_cut"], 94)
                ev = [line for line in (folder / "ssb_events.jsonl").read_text().splitlines() if line]
                self.assertEqual(len(ev), 6)
                i, q = pattern(600, "CAP16", 0)
                np.testing.assert_array_equal(iq[0].real, i)
        finally:
            dev.stop()


if __name__ == "__main__":
    unittest.main()
