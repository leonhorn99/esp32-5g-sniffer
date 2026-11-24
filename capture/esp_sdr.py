"""Client for the ESP-SDR serial protocol (ESPARGOS esp-sdr, protocol 6) and the
ESPsoup ESP32-C5 extensions.

Protocol facts used here (checked 2026-11-19 against esp-sdr docs/protocol.md and
espsoup-firmware docs/espsoup-c5.md, version 0.7.1):

* Line-based ASCII commands ending in ``\\n``; exactly one text reply line per
  command, replies starting with ``ERR`` are errors. Native USB or UART at 2 Mbaud.
* ``SYNC <nonce>`` is echoed; everything before the echo is discarded.
* ``CAP16|CAP20 <samples> <rate_code>`` -> ``DATA <n> <crc32-hex> <capture-us>``
  followed by ceil(n * bits_per_sample / 8) bytes (16 or 20 bits per sample).
* ESPsoup ``CAPT16|CAPT20`` -> ``DATAT <n> <crc32-hex> <us> <t0_us> <seq> <lag>``
  with the same payload; t0_us is the time of sample 0 in microseconds since boot.
* The firmware sends (I, -Q): Q is negated on the host (done in ``unpack_iq``).
* CAP16: per sample, byte 0 = I (int8), byte 1 = -Q (int8).
  CAP20: 20-bit samples packed LSB first without padding, I = bits 0-9,
  -Q = bits 10-19, both signed 10-bit.
* CRC32 is zlib's polynomial over the payload bytes.
"""
from __future__ import annotations

import json
import math
import secrets
import time
import zlib
from dataclasses import dataclass, field

import numpy as np

# Rate codes of CAP16/CAP20 (protocol.md). The C5 supports codes 0-5.
RATE_CODES = {
    0: 80e6, 1: 40e6, 2: 20e6, 3: 10e6, 4: 8e6, 5: 4e6,
    6: 16e6, 7: 32e6, 8: 32e6 / 3, 9: 6.4e6,
}
BITS_PER_SAMPLE = {"CAP16": 16, "CAP20": 20}


class EspSdrError(RuntimeError):
    """The firmware answered ERR, or the reply did not match the protocol."""


def rate_code_for(fs_hz: float) -> int:
    """Rate code for a nominal sample rate in Hz (e.g. 20e6 -> 2)."""
    for code, fs in RATE_CODES.items():
        if abs(fs - fs_hz) < 1.0:
            return code
    raise ValueError(f"no rate code for {fs_hz} Hz; known: {sorted(RATE_CODES.values())}")


def payload_bytes(n: int, fmt: str) -> int:
    return math.ceil(n * BITS_PER_SAMPLE[fmt] / 8)


def unpack_iq(payload: bytes, n: int, fmt: str) -> np.ndarray:
    """Decode a CAP16/CAP20 payload to complex64 (I + jQ) in ADC counts.

    Q is negated here, so the result is in standard orientation.
    """
    need = payload_bytes(n, fmt)
    if len(payload) < need:
        raise ValueError(f"payload has {len(payload)} bytes, {fmt} with n={n} needs {need}")
    buf = np.frombuffer(payload, dtype=np.uint8, count=need)
    if fmt == "CAP16":
        v = buf.view(np.int8).astype(np.float32).reshape(n, 2)
        i, mq = v[:, 0], v[:, 1]
    elif fmt == "CAP20":
        # Two samples = 40 bits = 5 bytes. Pad to whole groups of 5 bytes.
        groups = (n + 1) // 2
        b = np.zeros(groups * 5, dtype=np.uint64)
        b[:need] = buf
        b = b.reshape(groups, 5)
        word = b[:, 0] | (b[:, 1] << 8) | (b[:, 2] << 16) | (b[:, 3] << 24) | (b[:, 4] << 32)
        s = np.empty(groups * 2, dtype=np.uint64)
        s[0::2] = word & 0xFFFFF
        s[1::2] = (word >> 20) & 0xFFFFF
        s = s[:n]
        i = _sign10(s & 0x3FF)
        mq = _sign10((s >> 10) & 0x3FF)
    else:
        raise ValueError(f"unknown format {fmt!r}")
    return (i - 1j * mq).astype(np.complex64)


def raw_words(payload: bytes, n: int, fmt: str) -> np.ndarray:
    """The device's raw dump words (I = bits 0-9, -Q = bits 10-19) behind a payload.

    CAP20 carries bits 0-19 exactly. CAP16 carries bits 2-9 and 12-19, so the two
    lowest bits of each field come back as 0. Used to run the on-chip DSP on the host.
    """
    if fmt == "CAP20":
        groups = (n + 1) // 2
        b = np.zeros(groups * 5, dtype=np.uint64)
        need = payload_bytes(n, fmt)
        b[:need] = np.frombuffer(payload, dtype=np.uint8, count=need)
        b = b.reshape(groups, 5)
        word = b[:, 0] | (b[:, 1] << 8) | (b[:, 2] << 16) | (b[:, 3] << 24) | (b[:, 4] << 32)
        s = np.empty(groups * 2, dtype=np.uint64)
        s[0::2] = word & 0xFFFFF
        s[1::2] = (word >> 20) & 0xFFFFF
        return s[:n].astype(np.uint32)
    if fmt == "CAP16":
        v = np.frombuffer(payload, dtype=np.uint8, count=2 * n).astype(np.uint32).reshape(n, 2)
        return ((v[:, 0] << 2) & 0x3FF) | (((v[:, 1] << 2) & 0x3FF) << 10)
    raise ValueError(f"unknown format {fmt!r}")


def pack_iq(i: np.ndarray, q: np.ndarray, fmt: str) -> bytes:
    """Inverse of unpack_iq (device orientation: sends I and -Q). Used by tests."""
    i = np.asarray(i, dtype=np.int64)
    mq = -np.asarray(q, dtype=np.int64)
    n = len(i)
    if fmt == "CAP16":
        if i.min(initial=0) < -128 or i.max(initial=0) > 127 or mq.min(initial=0) < -128 or mq.max(initial=0) > 127:
            raise ValueError("CAP16 values must fit int8")
        out = np.empty(2 * n, dtype=np.int8)
        out[0::2], out[1::2] = i, mq
        return out.tobytes()
    if fmt == "CAP20":
        if i.min(initial=0) < -512 or i.max(initial=0) > 511 or mq.min(initial=0) < -512 or mq.max(initial=0) > 511:
            raise ValueError("CAP20 values must fit signed 10 bit")
        s = ((i & 0x3FF) | ((mq & 0x3FF) << 10)).astype(np.uint64)
        if n % 2:
            s = np.append(s, np.uint64(0))
        word = s[0::2] | (s[1::2] << 20)
        b = np.empty((len(word), 5), dtype=np.uint8)
        for k in range(5):
            b[:, k] = (word >> (8 * k)) & 0xFF
        return b.reshape(-1).tobytes()[: payload_bytes(n, fmt)]
    raise ValueError(f"unknown format {fmt!r}")


def _sign10(x: np.ndarray) -> np.ndarray:
    x = x.astype(np.int64)
    return ((x ^ 0x200) - 0x200).astype(np.float32)


@dataclass
class Burst:
    """One captured burst, with the raw payload kept as received."""
    fmt: str                 # "CAP16" or "CAP20"
    n: int
    rate_code: int
    payload: bytes
    crc_ok: bool
    crc32: int
    capture_us: int
    t0_us: int | None = None  # CAPT only
    seq: int | None = None
    lag: int | None = None
    reply: str = ""
    host_time: float = 0.0
    extra: dict = field(default_factory=dict)

    @property
    def fs(self) -> float:
        return RATE_CODES[self.rate_code]

    def iq(self) -> np.ndarray:
        return unpack_iq(self.payload, self.n, self.fmt)


@dataclass
class SsbFrame:
    """One SSBDET frame (firmware branch ssb-sniffer): the raw I/Q around one detected SSB.

    The first sample is at t0_us + t0_sub_ns / 1000 + (flip + start) / fs, sample j at raw
    index start + j R; pss_idx is the raw index (in the capture) of the PSS useful start, so
    (pss_idx - start) / R is its index in this cut. n = 0: header-only frame (the firmware
    sends the I/Q of every k-th SSB only). From firmware ssb3 on: R > 1 = decimated cut at
    fs / R with the SSB centre at 0 Hz and cfo_ddc_hz already removed (samples shifted right
    by gsh bits), pci from the on-chip SSS (-1 below threshold), sss its statistic.
    """
    seq: int
    n: int
    crc32: int
    t0_us: int
    flip: int
    start: int
    fs: int
    bits: int
    mode: str            # "A" acquisition burst, "T" tracking burst
    stat: float          # PSS statistic (noise: mean 1)
    nid2: int
    cfo_hz: int          # CFO measured on this SSB
    pss_idx: int
    payload: bytes
    crc_ok: bool
    host_time: float
    reply: str = ""
    t0_sub_ns: int = 0         # ssb3: sub-microsecond part of t0 (62.5 ns steps)
    pci: int | None = None     # ssb3: PCI from the on-chip SSS (-1: below threshold)
    sss: float | None = None   # ssb3: SSS statistic
    R: int = 1                 # ssb3: decimation of the cut (1 = raw)
    gsh: int = 0               # ssb3: right shift applied to a decimated cut
    cfo_ddc_hz: int = 0        # ssb3: CFO removed from a decimated cut

    @property
    def fmt(self) -> str:
        return "CAP16" if self.bits == 16 else "CAP20"

    @property
    def t_first_us(self) -> float:
        return self.t0_us + self.t0_sub_ns / 1000 + (self.flip + self.start) / self.fs * 1e6

    @property
    def fs_cut(self) -> float:
        """Sample rate of the payload."""
        return self.fs / self.R

    @property
    def pss_in_cut(self) -> float:
        return (self.pss_idx - self.start) / self.R

    def iq(self) -> np.ndarray:
        return unpack_iq(self.payload, self.n, self.fmt)


class EspSdr:
    """Talks to one board. ``port`` is a device path such as /dev/ttyACM0.

    Pass an already open serial-like object as ``ser`` for tests.
    """

    def __init__(self, port: str | None = None, baud: int = 2_000_000,
                 timeout: float = 3.0, ser=None):
        if ser is None:
            import serial  # pyserial
            # write_timeout: firmware that does not read the USB port (e.g. a factory demo)
            # would otherwise block the first write forever
            ser = serial.Serial(port, baudrate=baud, timeout=timeout, write_timeout=timeout)
        self.ser = ser
        self.timeout = timeout
        self._linebuf = bytearray()    # partial line kept across read timeouts
        self.caps: set[str] = set()
        self.info_line = ""
        self.limits: dict = {}

    # ----- low level -----------------------------------------------------
    def close(self):
        try:
            self.cmd("RELEASE")
        except Exception:
            pass
        try:
            self.ser.reset_output_buffer()
        except Exception:
            pass
        self.ser.close()

    def hard_reset(self, settle: float = 2.5):
        """Restart the chip through the native USB port (USB-Serial/JTAG reset logic, as esptool's
        hard reset): DTR low (normal boot, not download mode), RTS high = reset, RTS low = run.
        The USB connection stays up. Needed on 23 Nov: after resets of its own (USB re-enumeration)
        the firmware did not answer until reset this way."""
        self.ser.dtr = False
        self.ser.rts = True
        time.sleep(0.2)
        self.ser.rts = False
        time.sleep(settle)                      # bootloader + app start (boot log about 0.3 s)
        self.ser.reset_input_buffer()
        self._linebuf.clear()

    def _write_line(self, line: str):
        self.ser.write((line + "\n").encode("ascii"))
        self.ser.flush()

    def _read_line(self, timeout: float | None = None) -> str:
        """One text line. A timeout keeps the bytes received so far for the next call:
        dropping a half-received header (e.g. after a pause in a stream) would make the
        following binary payload be read as text (lost frames, SSBDET test of 23 Nov)."""
        deadline = time.monotonic() + (self.timeout if timeout is None else timeout)
        buf = self._linebuf
        while time.monotonic() < deadline:
            b = self.ser.read(1)
            if not b:
                continue
            if b == b"\n":
                line = buf.decode("ascii", errors="replace").rstrip("\r")
                buf.clear()
                return line
            buf += b
        raise TimeoutError(f"no reply line within timeout (partial so far {bytes(buf)[:80]!r})")

    def _read_exact(self, nbytes: int, timeout: float) -> bytes:
        deadline = time.monotonic() + timeout
        out = bytearray()
        while len(out) < nbytes and time.monotonic() < deadline:
            chunk = self.ser.read(nbytes - len(out))
            if chunk:
                out += chunk
        if len(out) < nbytes:
            raise TimeoutError(f"payload: got {len(out)} of {nbytes} bytes")
        return bytes(out)

    def cmd(self, line: str, timeout: float | None = None) -> str:
        """Send one command and return its reply line (raises on ERR)."""
        self._write_line(line)
        reply = self._read_line(timeout)
        if reply.startswith("ERR"):
            raise EspSdrError(f"{line!r} -> {reply!r}")
        return reply

    # ----- session -------------------------------------------------------
    def sync(self, attempts: int = 5) -> None:
        """Resynchronise: send SYNC <nonce> and discard everything up to the echo."""
        for _ in range(attempts):
            # decimal: the firmware parses "SYNC %" SCNu64 (receiver.c); a hex nonce gives "ERR command"
            nonce = str(secrets.randbelow(10**12))
            self.ser.reset_input_buffer()
            self._linebuf.clear()
            try:
                self._write_line(f"SYNC {nonce}")
            except Exception as e:  # serial.SerialTimeoutException: nobody reads on the board side
                self.ser.reset_output_buffer()   # else closing the tty waits 30 s for the stuck bytes
                raise EspSdrError(f"SYNC failed: the board does not read its USB serial port ({e}); "
                                  "it is not running ESP-SDR/ESPsoup firmware (flash it: firmware/fw.sh flash)")
            deadline = time.monotonic() + self.timeout
            while time.monotonic() < deadline:
                try:
                    line = self._read_line(timeout=max(0.05, deadline - time.monotonic()))
                except TimeoutError:
                    break
                if line.strip().endswith(f"SYNC {nonce}"):
                    return
        raise EspSdrError("SYNC failed: board not answering (wrong port, board resetting, or ModemManager?)")

    def connect(self) -> dict:
        """SYNC, then read INFO, CAPS, LIMITS?, VERSION? and FWOPT? (if offered)."""
        self.sync()
        self.info_line = self.cmd("INFO")
        caps = self.cmd("CAPS")
        self.caps = set(caps.split()[1:]) if caps.startswith("CAPS") else set()
        out = {"info": self.info_line, "caps": sorted(self.caps)}
        try:  # protocol 6 advertises RXLIMITS; older firmware answers ERR
            out["limits"] = self.limits = _json_after_keyword(self.cmd("LIMITS?"), "LIMITS")
        except (EspSdrError, ValueError):
            out["limits"] = None
        for q, key in (("VERSION?", "VERSION"),):
            try:
                out["version"] = _json_after_keyword(self.cmd(q), key)
            except (EspSdrError, ValueError):
                out["version"] = None
        if "FWOPT" in self.caps:
            out["fwopt"] = self.cmd("FWOPT?")
        if "BOOT" in self.caps:            # ssb3: why the chip last started (brownout, panic, ...)
            out["boot"] = self._try("BOOT?")
        out["transport"] = self._try("TRANSPORT?")
        out["range"] = self._try("RANGE?")
        return out

    def _try(self, line: str) -> str | None:
        try:
            return self.cmd(line)
        except (EspSdrError, TimeoutError):
            return None

    @property
    def max_samples(self) -> int:
        if self.limits and "maxSamples" in self.limits:
            return int(self.limits["maxSamples"])
        parts = self.info_line.split()
        return int(parts[-1]) if parts and parts[-1].isdigit() else 16380

    @property
    def family(self) -> str:
        return self.info_line.split()[0] if self.info_line else ""

    # ----- receiver settings ---------------------------------------------
    def set_freq_mhz(self, mhz: int) -> str:
        if int(mhz) != mhz:
            raise ValueError("FREQ takes whole MHz; shift the signal digitally on the host")
        return self.cmd(f"FREQ {int(mhz)}")

    def set_gain(self, index: int | None) -> str:
        """Manual gain index (PHY table entry, not dB); None = hardware AGC."""
        return self.cmd("GAIN HARDWARE" if index is None else f"GAIN MANUAL {int(index)}")

    def get_gain(self) -> str:
        return self.cmd("GAIN?")

    def set_bandwidth_mhz(self, mhz: int) -> str:
        """Approximate analog bandwidth; 0 = widest (C5: 11-48 MHz)."""
        return self.cmd(f"BANDWIDTH {int(mhz)}")

    def temperature(self) -> float | None:
        if "TEMP" not in self.caps:
            return None
        r = self._try("TEMP?")
        try:
            return float(r.split()[1]) if r else None
        except (IndexError, ValueError):
            return None

    def dc_calibrate(self) -> str | None:
        """ESPsoup hardware DC calibration (needs manual gain; lasts until the next FREQ/GAIN)."""
        if "DCCAL" not in self.caps:
            return None
        return self.cmd("DCCAL", timeout=self.timeout + 2.0)

    # ----- capture -------------------------------------------------------
    def capture(self, n: int, rate_code: int, bits: int = 10, timed: bool = True) -> Burst:
        """One burst. bits=8 -> CAP16, bits=10 -> CAP20. timed uses CAPT if offered."""
        fmt = {8: "CAP16", 10: "CAP20"}[bits]
        use_t = timed and "CAPT" in self.caps
        command = f"{'CAPT' if use_t else 'CAP'}{fmt[3:]} {int(n)} {int(rate_code)}"
        self._write_line(command)
        host_time = time.time()
        reply = self._read_line()
        if reply.startswith("ERR"):
            raise EspSdrError(f"{command!r} -> {reply!r}")
        parts = reply.split()
        if use_t:
            if parts[0] != "DATAT" or len(parts) < 7:
                raise EspSdrError(f"unexpected reply to {command!r}: {reply!r}")
            n_rx, crc, us, t0, seq, lag = int(parts[1]), int(parts[2], 16), int(parts[3]), int(parts[4]), int(parts[5]), int(parts[6])
        else:
            if parts[0] != "DATA" or len(parts) < 4:
                raise EspSdrError(f"unexpected reply to {command!r}: {reply!r}")
            n_rx, crc, us = int(parts[1]), int(parts[2], 16), int(parts[3])
            t0 = seq = lag = None
        nbytes = payload_bytes(n_rx, fmt)
        # USB full speed carries about 0.93 MB/s; allow 4x that time plus the timeout.
        payload = self._read_exact(nbytes, timeout=self.timeout + 4 * nbytes / 0.9e6)
        crc_ok = (zlib.crc32(payload) & 0xFFFFFFFF) == crc
        return Burst(fmt=fmt, n=n_rx, rate_code=rate_code, payload=payload, crc_ok=crc_ok,
                     crc32=crc, capture_us=us, t0_us=t0, seq=seq, lag=lag, reply=reply,
                     host_time=host_time)


def ssb_stream(sdr: "EspSdr", rate_code: int, offset_hz: int, nid2: int, cfo_hz: int, thr10: int = 140,
               pre_us: int = 20, post_us: int = 20, period_us: int = 20000, bits: int = 8, every: int = 1,
               max_frames: int = 0, track: bool = True, duration: float | None = None, stop=None,
               idle_timeout: float | None = None, dec: bool = False, cfo_span: int = 0):
    """Run SSBDET (on-chip SSB detection) and yield ("start", line), ("frame", SsbFrame)...,
    ("end", line). Other text lines from the board (ESP-IDF log, e.g. a brownout or panic
    message before a reset) come as ("text", line). The run ends after max_frames, after
    duration seconds, or when stop() returns True (both send a stop line and wait for the END line).
    Firmware ssb3 (CAPS SSBDET3): nid2 3 searches all three N_ID2, cfo_span > 0 searches the CFO
    +- cfo_span in 5 kHz steps, dec sends cuts decimated to 5 MS/s (about 4x less data).
    idle_timeout: after that many seconds without any line, send STOP: an END line then comes as
    ("idle", line) (alive, but no SSB found), no END within 5 s raises TimeoutError (board dead)."""
    fmt_bits = 16 if bits == 8 else 20
    command = (f"SSBDET {rate_code} {int(offset_hz)} {nid2} {int(cfo_hz)} {thr10} {pre_us} {post_us} "
               f"{period_us} {fmt_bits} {every} {max_frames} {1 if track else 0}")
    if dec or cfo_span or nid2 > 2:
        if "SSBDET3" not in sdr.caps:
            raise EspSdrError("cell search and decimated cuts need firmware ssb3 (CAPS SSBDET3): flash the current build")
        command += f" {1 if dec else 0} {int(cfo_span)}"
    sdr._write_line(command)
    reply = sdr._read_line(timeout=sdr.timeout + 5)
    if not reply.startswith("OK SSBDET"):
        raise EspSdrError(f"{command!r} -> {reply!r}")
    yield "start", reply
    t_end = time.monotonic() + duration if duration else None
    stop_sent, idle_stop, t_stop = False, False, 0.0
    last_rx = time.monotonic()
    while True:
        if not stop_sent and ((t_end and time.monotonic() > t_end) or (stop and stop())):
            sdr._write_line("STOP")
            stop_sent, t_stop = True, time.monotonic()
        if idle_timeout and not stop_sent and time.monotonic() - last_rx > idle_timeout:
            sdr._write_line("STOP")
            stop_sent, idle_stop, t_stop = True, True, time.monotonic()
        if stop_sent and time.monotonic() - t_stop > 5.0:
            raise TimeoutError("SSBDET: no END line within 5 s after STOP (board reset or hung?)")
        try:
            line = sdr._read_line(timeout=1.0)
        except TimeoutError:
            continue
        last_rx = time.monotonic()
        if line.startswith("END SSBDET"):
            yield ("idle" if idle_stop else "end"), line
            return
        if line.startswith("ERR"):
            raise EspSdrError(f"SSBDET: {line!r}")
        if not line.startswith("SSBF "):
            if line.strip():
                yield "text", line
            continue
        p = line.split()
        # SSBF seq n crc t0_us flip start fs bits mode stat_x10 nid2 cfo_hz pss_idx
        #      [t0_sub_ns pci sss_x10 R gsh cfo_ddc_hz]   (ssb3)
        n, b = int(p[2]), int(p[8])
        ext = dict(t0_sub_ns=int(p[14]), pci=int(p[15]), sss=int(p[16]) / 10, R=int(p[17]), gsh=int(p[18]),
                   cfo_ddc_hz=int(p[19])) if len(p) >= 20 else {}
        fmt = "CAP16" if b == 16 else "CAP20"
        host_time = time.time()
        payload = sdr._read_exact(payload_bytes(n, fmt), timeout=sdr.timeout + 4 * payload_bytes(n, fmt) / 0.9e6) if n else b""
        crc = int(p[3], 16)
        yield "frame", SsbFrame(seq=int(p[1]), n=n, crc32=crc, t0_us=int(p[4]), flip=int(p[5]), start=int(p[6]),
                                fs=int(p[7]), bits=b, mode=p[9], stat=int(p[10]) / 10, nid2=int(p[11]),
                                cfo_hz=int(p[12]), pss_idx=int(p[13]), payload=payload,
                                crc_ok=(not n) or (zlib.crc32(payload) & 0xFFFFFFFF) == crc,
                                host_time=host_time, reply=line, **ext)


def _json_after_keyword(reply: str, keyword: str):
    if not reply.startswith(keyword):
        raise ValueError(f"expected {keyword}, got {reply!r}")
    return json.loads(reply[len(keyword):].strip())
