"""A fake ESP-SDR/ESPsoup C5 board on a pseudo-terminal, for tests without hardware.

It answers the commands the capture client uses and sends a deterministic sample
pattern, so decoders can be checked value by value:
    I[k] = ((37*k + 11*seq) mod M) - M/2,  Q[k] = ((91*k + 5*seq) mod (M-1)) - (M/2 - 1)
with M = 256 (CAP16) or 1024 (CAP20); Q stops at -(M/2-1) so that -Q fits.
The device sends (I, -Q) like the real firmware.
"""
from __future__ import annotations

import os
import threading
import tty
import zlib

import numpy as np

from capture.esp_sdr import pack_iq


def pattern(n: int, fmt: str, seq: int = 0):
    m = 256 if fmt == "CAP16" else 1024
    k = np.arange(n)
    i = (37 * k + 11 * seq) % m - m // 2
    q = (91 * k + 5 * seq) % (m - 1) - (m // 2 - 1)
    return i, q


class FakeDevice(threading.Thread):
    def __init__(self, caps=("RXLIMITS", "TUNEEXT", "GPIO", "BIGCAP", "CAPT", "TEMP", "FWOPT", "DCCAL"),
                 max_samples=32704, junk_before_sync=True, corrupt_seqs=(), gain_scaling=False):
        super().__init__(daemon=True)
        self.master, self.slave = os.openpty()
        tty.setraw(self.slave)
        self.port = os.ttyname(self.slave)
        self.caps = list(caps)
        self.max_samples = max_samples
        self.junk = junk_before_sync
        self.corrupt_seqs = set(corrupt_seqs)
        self.gain_scaling = gain_scaling   # scale the pattern by gain/80 (for auto-gain tests)
        self.gain = 80
        self.seq = 0
        self.t_us = 1_000_000
        self.log: list[str] = []
        self._stop = threading.Event()

    def _send(self, data: bytes):
        os.write(self.master, data)

    def _line(self, s: str):
        self._send((s + "\n").encode())

    def run(self):
        buf = b""
        if self.junk:  # a board that just booted prints a log line and some binary
            self._send(b"rst:0x15 boot log\n\x00\xff\x13")
        while not self._stop.is_set():
            try:
                data = os.read(self.master, 4096)
            except OSError:
                return
            if not data:
                continue
            buf += data
            while b"\n" in buf:
                line, buf = buf.split(b"\n", 1)
                self.handle(line.decode().strip())

    def stop(self):
        self._stop.set()
        for fd in (self.master, self.slave):
            try:
                os.close(fd)
            except OSError:
                pass

    def _ssb_emit(self, fs, bits, max_frames):
        """SSBDET run: one frame every 5 ms (pattern payload of 600 samples) until stopped."""
        fmt = "CAP16" if bits == 16 else "CAP20"
        k = 0
        while not self._ssb_stop.is_set() and (not max_frames or k < max_frames):
            n = 0 if k % 3 == 2 else 600              # every third frame header-only
            i, q = pattern(n, fmt, k)
            payload = pack_iq(i, q, fmt) if n else b""
            crc = zlib.crc32(payload) & 0xFFFFFFFF if n else 0
            mode = "A" if k == 0 else "T"
            self._line(f"SSBF {k + 1} {n} {crc:08x} {self.t_us} 0 1000 {fs} {bits} {mode} 215 0 11250 1094")
            if n:
                self._send(payload)
            self.t_us += 20_000
            k += 1
            self._ssb_stop.wait(0.005)
        self._line(f"END SSBDET {k} {k} 1 {max(k - 1, 0)} 0 60 900 4000 300 11250 {'limit' if max_frames and k >= max_frames else 'stop'}")
        self._ssb_active = False

    def handle(self, line: str):
        self.log.append(line)
        if not line:
            return
        if getattr(self, "_ssb_active", False):    # any line stops a running SSBDET
            self._ssb_stop.set()
            return
        w = line.split()
        c = w[0]
        if c == "SYNC":   # like the firmware: sscanf(line, "SYNC %" SCNu64 " %c") == 1
            if len(w) == 2 and w[1].isdigit():
                self._line(f"SYNC {int(w[1])}")
            else:
                self._line("ERR command")
        elif c == "INFO":
            self._line("C5SDR 6 burst 16380")
        elif c == "CAPS":
            self._line("CAPS " + " ".join(self.caps))
        elif c == "LIMITS?":
            self._line('LIMITS {"gain":[0,80,1],"bandwidth":[11,48,1,48],"rates":[80000000,40000000,'
                       '20000000,10000000,8000000,4000000],"bits":[8,10],"maxSamples":%d}' % self.max_samples)
        elif c == "VERSION?":
            self._line('VERSION {"git":"fake","built":"2026-11-19"}')
        elif c == "FWOPT?" and "FWOPT" in self.caps:
            self._line("FWOPT 0.7.1-dev espsoup rx-only C5 (fake)")
        elif c == "TRANSPORT?":
            self._line("TRANSPORT USB 0")
        elif c == "RANGE?":
            self._line("RANGE 100 6000 1")
        elif c in ("RELEASE", "FREQ", "BANDWIDTH"):
            self._line("OK")
        elif c == "GAIN":
            if len(w) == 3 and w[1] == "MANUAL":
                self.gain = int(w[2])
            self._line("OK")
        elif c == "GAIN?":
            self._line("GAIN MANUAL 50")
        elif c == "TEMP?" and "TEMP" in self.caps:
            self._line("TEMP 43.9")
        elif c == "DCCAL" and "DCCAL" in self.caps:
            self._line("DCCAL before 2.1 -1.3 after 0.3 0.2 base 0 0 best 0 0 steps 3 jac 1 0 0 1 auto 0")
        elif c == "SSBDET" and "SSBDET" in self.caps:
            a = [int(x) for x in w[1:]]
            rates = {0: 80_000_000, 1: 40_000_000, 2: 20_000_000, 4: 8_000_000, 5: 4_000_000}
            if len(a) < 4 or a[0] not in rates:
                self._line("ERR ssbdet_args")
                return
            fs = rates[a[0]]
            bits = a[8] if len(a) > 8 else 16
            max_frames = a[10] if len(a) > 10 else 0
            self._line(f"OK SSBDET {fs} {fs // 4_000_000} {6 * (fs // 4_000_000) + 1} 493 140 20000")
            self._ssb_stop = threading.Event()
            self._ssb_active = True
            threading.Thread(target=self._ssb_emit, args=(fs, bits, max_frames), daemon=True).start()
        elif c in ("CAP16", "CAP20", "CAPT16", "CAPT20"):
            timed = c.startswith("CAPT")
            if timed and "CAPT" not in self.caps:
                self._line("ERR unknown_command")
                return
            fmt = "CAP" + c[-2:]
            n, rate = int(w[1]), int(w[2])
            if n > self.max_samples or rate not in range(6):
                self._line("ERR cap_args")
                return
            i, q = pattern(n, fmt, self.seq)
            if self.gain_scaling:
                i, q = (np.round(i * self.gain / 80).astype(int), np.round(q * self.gain / 80).astype(int))
            payload = pack_iq(i, q, fmt)
            crc = zlib.crc32(payload) & 0xFFFFFFFF
            if self.seq in self.corrupt_seqs:
                payload = bytes([payload[0] ^ 0x01]) + payload[1:]
            us = 821
            if timed:
                self._line(f"DATAT {n} {crc:08x} {us} {self.t_us} {self.seq} 7")
            else:
                self._line(f"DATA {n} {crc:08x} {us}")
            self._send(payload)
            self.seq += 1
            self.t_us += 80_000
        else:
            self._line("ERR unknown_command")


def write_fixture(folder):
    """Write the small session used by tests/test_matlab_loader.m (CAP20 odd n, CAP16)."""
    from capture.esp_sdr import EspSdr
    from capture.session import SessionWriter

    dev = FakeDevice()
    dev.start()
    try:
        sdr = EspSdr(dev.port, timeout=2.0)
        device = sdr.connect()
        w = SessionWriter(folder, device, {"note": "fixture from tests/fake_device.py, pattern() values"})
        b = sdr.capture(101, 2, bits=10)          # seq 0, CAPT20, odd length
        w.add(b, freq_mhz=2138, gain=50, bandwidth_mhz=0)
        b = sdr.capture(64, 5, bits=8)            # seq 1, CAPT16
        w.add(b, freq_mhz=2138, gain=50, bandwidth_mhz=0)
        w.close()
        sdr.ser.close()
    finally:
        dev.stop()


if __name__ == "__main__":
    import sys
    write_fixture(sys.argv[1])
