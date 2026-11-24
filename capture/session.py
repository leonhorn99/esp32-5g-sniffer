"""Session files: the raw burst payloads plus one JSON line of metadata per burst.

A session is a folder with
  session.json   device answers (INFO, CAPS, LIMITS, VERSION, FWOPT) and the settings
  bursts.bin     all payloads exactly as received, appended
  bursts.jsonl   one line per burst: byte offset/length in bursts.bin, format, n,
                 rate, frequency, gain, timestamps, CRC result, temperature

The MATLAB reader is matlab/esp_load_session.m; both readers unpack the same way.
"""
from __future__ import annotations

import json
import os
import platform
import time
from pathlib import Path

import numpy as np

from .esp_sdr import RATE_CODES, Burst, unpack_iq

SESSION_FORMAT = 1


class SessionWriter:
    def __init__(self, folder: str | Path, device: dict, settings: dict):
        self.folder = Path(folder)
        self.folder.mkdir(parents=True, exist_ok=False)
        meta = {
            "format": SESSION_FORMAT,
            "created_unix": time.time(),
            "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "host": platform.node(),
            "device": device,
            "settings": settings,
        }
        (self.folder / "session.json").write_text(json.dumps(meta, indent=2))
        self._bin = open(self.folder / "bursts.bin", "ab")
        self._jsonl = open(self.folder / "bursts.jsonl", "a")
        self._offset = self._bin.tell()
        self.count = 0

    def add(self, burst: Burst, freq_mhz: float, gain, bandwidth_mhz, temp_c=None, extra=None):
        self._bin.write(burst.payload)
        rec = {
            "idx": self.count,
            "offset": self._offset,
            "nbytes": len(burst.payload),
            "fmt": burst.fmt,
            "n": burst.n,
            "rate_code": burst.rate_code,
            "fs": burst.fs,
            "freq_mhz": freq_mhz,
            "gain": gain,
            "bandwidth_mhz": bandwidth_mhz,
            "crc_ok": burst.crc_ok,
            "crc32": f"{burst.crc32:08x}",
            "capture_us": burst.capture_us,
            "t0_us": burst.t0_us,
            "seq": burst.seq,
            "lag": burst.lag,
            "host_time": burst.host_time,
            "temp_c": temp_c,
            "reply": burst.reply,
        }
        if extra:
            rec.update(extra)
        self._jsonl.write(json.dumps(rec) + "\n")
        self._offset += len(burst.payload)
        self.count += 1
        if self.count % 20 == 0:
            self.flush()
        return rec

    def flush(self):
        self._bin.flush()
        self._jsonl.flush()
        os.fsync(self._bin.fileno())

    def close(self):
        self.flush()
        self._bin.close()
        self._jsonl.close()


class SsbRecorder:
    """Session of SSBDET frames (on-chip SSB detection): every cut with raw I/Q becomes a burst
    (t0_us = time of the cut's first sample, plus the ssb_* fields), every frame, also the
    header-only ones, a line of ssb_events.jsonl. Used by capture_bursts ssb-stream and ssb_live."""

    def __init__(self, folder: str | Path, device: dict, settings: dict, freq_mhz: float, gain,
                 bandwidth_mhz, rate_code: int):
        self.writer = SessionWriter(folder, device, settings)
        self.folder = self.writer.folder
        self.events = open(self.folder / "ssb_events.jsonl", "a")
        self.freq_mhz, self.gain, self.bandwidth_mhz, self.rate_code = freq_mhz, gain, bandwidth_mhz, rate_code

    def add(self, f) -> dict:
        """f: esp_sdr.SsbFrame. Returns the event line written."""
        ev = {"seq": f.seq, "n": f.n, "t0_us": f.t0_us, "flip": f.flip, "start": f.start, "fs": f.fs,
              "mode": f.mode, "stat": f.stat, "nid2": f.nid2, "cfo_hz": f.cfo_hz, "pss_idx": f.pss_idx,
              "t_pss_us": f.t0_us + (f.flip + f.pss_idx) / f.fs * 1e6, "crc_ok": f.crc_ok, "host_time": f.host_time}
        if f.n:
            b = Burst(fmt=f.fmt, n=f.n, rate_code=self.rate_code, payload=f.payload, crc_ok=f.crc_ok, crc32=f.crc32,
                      capture_us=0, t0_us=round(f.t_first_us, 3), seq=f.seq, lag=f.flip, reply=f.reply,
                      host_time=f.host_time)
            rec = self.writer.add(b, freq_mhz=self.freq_mhz, gain=self.gain, bandwidth_mhz=self.bandwidth_mhz,
                                  extra={"ssb_mode": f.mode, "ssb_stat": f.stat, "ssb_nid2": f.nid2,
                                         "ssb_cfo_hz": f.cfo_hz, "ssb_pss_in_cut": f.pss_idx - f.start,
                                         "capture_t0_us": f.t0_us, "cut_start": f.start})
            ev["burst_idx"] = rec["idx"]
        self.events.write(json.dumps(ev) + "\n")
        return ev

    @property
    def count(self) -> int:
        return self.writer.count

    def close(self):
        self.events.close()
        self.writer.close()


def load_session(folder: str | Path, only_crc_ok: bool = True):
    """Return (meta, records, iq_list). iq in ADC counts, Q already negated."""
    folder = Path(folder)
    meta = json.loads((folder / "session.json").read_text())
    records = [json.loads(line) for line in (folder / "bursts.jsonl").read_text().splitlines() if line.strip()]
    if only_crc_ok:
        records = [r for r in records if r["crc_ok"]]
    iq = []
    with open(folder / "bursts.bin", "rb") as f:
        for r in records:
            f.seek(r["offset"])
            iq.append(unpack_iq(f.read(r["nbytes"]), r["n"], r["fmt"]))
    return meta, records, iq


def fs_of(rate_code: int) -> float:
    return RATE_CODES[rate_code]


def burst_stats(iq: np.ndarray, fmt: str) -> dict:
    """RMS, peak and the share of samples at full scale (clipping indicator)."""
    full = 127 if fmt == "CAP16" else 511
    a = np.concatenate([iq.real, iq.imag])
    return {
        "rms": float(np.sqrt(np.mean(np.abs(iq) ** 2))),
        "peak": float(np.max(np.abs(a))),
        "clip_frac": float(np.mean(np.abs(a) >= full)),
        "dc": complex(np.mean(iq)),
    }
