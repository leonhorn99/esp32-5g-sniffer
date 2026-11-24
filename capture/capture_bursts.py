"""Capture I/Q bursts from an ESP32 running ESP-SDR / ESPsoup firmware.

Run from the project root:

  .venv_local/bin/python -m capture.capture_bursts info --port /dev/ttyACM0
  .venv_local/bin/python -m capture.capture_bursts gain-scan --port /dev/ttyACM0 --freq 2138 --rate 20
  .venv_local/bin/python -m capture.capture_bursts capture --port /dev/ttyACM0 \
        --freq 2138 --rate 20 --bits 10 --gain 50 --count 300 --label cable_m70dBm

Captures go to data/<YYYYmmdd_HHMMSS>_<label>/ (see capture/session.py for the format).
Use manual gain: the Wi-Fi AGC settles on packet preambles and clips on continuous signals.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from .esp_sdr import EspSdr, EspSdrError, rate_code_for
from .session import SessionWriter, SsbRecorder, burst_stats

ROOT = Path(__file__).resolve().parents[1]


def _open(args) -> tuple[EspSdr, dict]:
    sdr = EspSdr(args.port, baud=args.baud, timeout=args.timeout)
    device = sdr.connect()
    return sdr, device


def cmd_info(args):
    sdr, device = _open(args)
    try:
        device["temp_c"] = sdr.temperature()
        device["gain"] = sdr._try("GAIN?")
        print(json.dumps(device, indent=2))
    finally:
        sdr.close()


def _setup(sdr: EspSdr, freq_mhz: int, args, gain):
    sdr.set_freq_mhz(freq_mhz)
    if args.bandwidth is not None:
        sdr.set_bandwidth_mhz(args.bandwidth)
    sdr.set_gain(gain)
    if args.dccal and gain is not None:
        r = sdr.dc_calibrate()
        if r:
            print(f"  DCCAL: {r}")


def cmd_gain_scan(args):
    sdr, device = _open(args)
    try:
        g = (device.get("limits") or {}).get("gain") or [0, 80, 1]
        lo, hi = g[0], g[1]
        gains = range(lo, hi + 1, args.step) if args.gains is None else args.gains
        rate = rate_code_for(args.rate * 1e6)
        print(f"{sdr.family}  gain range {g}  FREQ {args.freq[0]} MHz  rate {args.rate} MS/s")
        print(" gain     rms    peak  clip%   |dc|")
        for gain in gains:
            _setup(sdr, args.freq[0], args, gain)
            b = sdr.capture(min(args.n, sdr.max_samples), rate, bits=args.bits)
            if not b.crc_ok:
                print(f"{gain:5d}  CRC error")
                continue
            s = burst_stats(b.iq(), b.fmt)
            print(f"{gain:5d} {s['rms']:7.2f} {s['peak']:7.1f} {100*s['clip_frac']:6.2f} {abs(s['dc']):6.2f}")
        print("Pick the highest gain with clip% = 0 and rms well below full scale "
              f"({127 if args.bits == 8 else 511}); with a strong cable signal aim for peak at about half scale.")
    finally:
        sdr.close()


def pick_gain(sdr: EspSdr, freq_mhz: int, rate: int, bits: int, args, step: int = 4,
              headroom: float = 0.5) -> tuple[int, list]:
    """Highest manual gain index with no clipped sample and peak <= headroom * full scale.

    One 8192-sample burst per gain index; returns (gain, scan table). If even the lowest
    gain exceeds the limit, the lowest gain is returned.
    """
    g = (sdr.limits or {}).get("gain") or [0, 80, 1]
    full = 127 if bits == 8 else 511
    table = []
    for gain in range(g[0], g[1] + 1, step):
        _setup(sdr, freq_mhz, args, gain)
        b = sdr.capture(min(8192, sdr.max_samples), rate, bits=bits)
        if not b.crc_ok:
            continue
        s = burst_stats(b.iq(), b.fmt)
        table.append({"gain": gain, "rms": s["rms"], "peak": s["peak"], "clip_frac": s["clip_frac"]})
    good = [r["gain"] for r in table if r["clip_frac"] == 0 and r["peak"] <= headroom * full]
    chosen = max(good) if good else (table[0]["gain"] if table else g[0])
    return chosen, table


def cmd_capture(args):
    sdr, device = _open(args)
    rate = rate_code_for(args.rate * 1e6)
    n = min(args.n or sdr.max_samples, sdr.max_samples)
    gain_scan = None
    if args.gain == "auto":
        args.gain, gain_scan = pick_gain(sdr, args.freq[0], rate, args.bits, args)
        print(f"auto gain: {args.gain} (highest index without clipping, peak <= half scale)")
    elif args.gain is not None:
        args.gain = int(args.gain)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    folder = Path(args.out) / (f"{stamp}_{args.label}" if args.label else stamp)
    settings = {
        "freq_mhz": args.freq, "per_freq": args.per_freq, "rate_msps": args.rate, "rate_code": rate,
        "bits": args.bits, "n": n, "gain": args.gain, "bandwidth_mhz": args.bandwidth,
        "dccal": args.dccal, "count": args.count, "duration_s": args.duration, "label": args.label,
        "note": args.note, "timed": "CAPT" in sdr.caps, "gain_scan": gain_scan,
    }
    writer = SessionWriter(folder, device, settings)
    print(f"{sdr.family} -> {folder}  ({n} samples, {args.rate} MS/s = {1e3*n/(args.rate*1e6):.2f} ms per burst)")
    t_start = time.monotonic()
    n_bad = 0
    k = 0
    fi = 0
    try:
        _setup(sdr, args.freq[0], args, args.gain)
        while True:
            if args.count and writer.count >= args.count:
                break
            if args.duration and time.monotonic() - t_start >= args.duration:
                break
            if len(args.freq) > 1 and k > 0 and k % args.per_freq == 0:
                fi = (fi + 1) % len(args.freq)
                _setup(sdr, args.freq[fi], args, args.gain)
            try:
                b = sdr.capture(n, rate, bits=args.bits)
            except (TimeoutError, EspSdrError) as e:
                print(f"  burst {k}: {e}; resync")
                sdr.sync()
                n_bad += 1
                k += 1
                continue
            temp = sdr.temperature() if (args.temp_every and k % args.temp_every == 0) else None
            st = burst_stats(b.iq(), b.fmt) if b.crc_ok else {}
            writer.add(b, freq_mhz=args.freq[fi], gain=args.gain, bandwidth_mhz=args.bandwidth,
                       temp_c=temp, extra={"rms": st.get("rms"), "clip_frac": st.get("clip_frac")})
            n_bad += 0 if b.crc_ok else 1
            k += 1
            if k % 25 == 0:
                el = time.monotonic() - t_start
                print(f"  {writer.count} bursts, {writer.count/el:.1f}/s, CRC/timeout errors {n_bad}, "
                      f"rms {st.get('rms', float('nan')):.1f}, clip {100*st.get('clip_frac', 0):.2f} %")
    except KeyboardInterrupt:
        print("stopped")
    finally:
        writer.close()
        sdr.close()
    print(f"done: {writer.count} bursts in {time.monotonic()-t_start:.1f} s, {n_bad} errors -> {folder}")


def cmd_ssb_stream(args):
    """On-chip SSB detection (firmware branch ssb-sniffer, SSBDET): store only the raw I/Q
    around each detected SSB. Cuts become bursts of a normal session (t0_us = time of the
    cut's first sample), every frame also goes to ssb_events.jsonl."""
    from .esp_sdr import ssb_stream
    sdr, device = _open(args)
    if "SSBDET" not in sdr.caps:
        sdr.close()
        raise SystemExit("this firmware has no SSBDET (flash the ssb-sniffer build: firmware/fw.sh flash)")
    rate = rate_code_for(args.rate * 1e6)
    gain_scan = None
    if args.gain == "auto":
        args.gain, gain_scan = pick_gain(sdr, args.freq[0], rate, 10, args)
        print(f"auto gain: {args.gain}")
    elif args.gain is not None:
        args.gain = int(args.gain)
    offset = int(round((args.ssb - args.freq[0]) * 1e6))
    stamp = time.strftime("%Y%m%d_%H%M%S")
    folder = Path(args.out) / (f"{stamp}_{args.label}" if args.label else f"{stamp}_ssbdet")
    settings = {"mode": "ssbdet", "freq_mhz": args.freq[:1], "ssb_mhz": args.ssb, "offset_hz": offset,
                "nid2": args.nid2, "cfo_hz": args.cfo, "rate_msps": args.rate, "rate_code": rate, "bits": args.bits,
                "gain": args.gain, "bandwidth_mhz": args.bandwidth, "thr": args.thr, "pre_us": args.pre,
                "post_us": args.post, "period_ms": args.period, "every": args.every, "track": not args.no_track,
                "count": args.count, "duration_s": args.duration, "label": args.label, "note": args.note,
                "gain_scan": gain_scan}
    rec = SsbRecorder(folder, device, settings, args.freq[0], args.gain, args.bandwidth, rate)
    _setup(sdr, args.freq[0], args, args.gain)
    t_start = time.monotonic()
    nframes, modes, bad = 0, {"A": 0, "T": 0}, 0
    stop = {"now": False}
    try:
        for kind, item in ssb_stream(sdr, rate, offset, args.nid2, args.cfo, int(round(args.thr * 10)), args.pre,
                                     args.post, int(args.period * 1000), args.bits, args.every, args.count,
                                     not args.no_track, args.duration, stop=lambda: stop["now"]):
            if kind == "start":
                print(item)
                continue
            if kind == "end":
                print(item)
                break
            if kind == "text":
                print(f"  board: {item}")
                continue
            f = item
            nframes += 1
            modes[f.mode] = modes.get(f.mode, 0) + 1
            bad += 0 if f.crc_ok else 1
            rec.add(f)
            if nframes % 25 == 0:
                el = time.monotonic() - t_start
                print(f"  {nframes} SSBs in {el:.1f} s ({nframes / el:.1f}/s), acquisition {modes['A']}, "
                      f"tracking {modes['T']}, CRC errors {bad}, last statistic {f.stat:.1f}, CFO {f.cfo_hz} Hz")
    except KeyboardInterrupt:
        print("stopping")
    finally:
        rec.close()
        sdr.close()
    el = time.monotonic() - t_start
    print(f"done: {nframes} SSB frames in {el:.1f} s, {rec.count} cuts stored -> {folder}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode", required=True)

    def common(sp):
        sp.add_argument("--port", required=True, help="e.g. /dev/ttyACM0 (native USB) or /dev/ttyUSB0")
        sp.add_argument("--baud", type=int, default=2_000_000)
        sp.add_argument("--timeout", type=float, default=3.0)

    def rx(sp):
        sp.add_argument("--freq", type=int, nargs="+", required=True, help="LO in whole MHz; several = round robin")
        sp.add_argument("--rate", type=float, default=20, help="MS/s: 80, 40, 20, 10, 8 or 4 on the C5")
        sp.add_argument("--bits", type=int, choices=(8, 10), default=10)
        sp.add_argument("--bandwidth", type=int, default=None, help="analog bandwidth in MHz, 0 = widest")
        sp.add_argument("--dccal", action="store_true", help="ESPsoup DCCAL after each retune (manual gain only)")

    s = sub.add_parser("info", help="print INFO/CAPS/LIMITS/VERSION")
    common(s)
    s.set_defaults(func=cmd_info)

    s = sub.add_parser("gain-scan", help="one burst per gain index: rms, peak, clipping")
    common(s)
    rx(s)
    s.add_argument("--n", type=int, default=16380)
    s.add_argument("--step", type=int, default=4)
    s.add_argument("--gains", type=int, nargs="*", default=None)
    s.set_defaults(func=cmd_gain_scan)

    s = sub.add_parser("capture", help="capture bursts into a session folder")
    common(s)
    rx(s)
    s.add_argument("--gain", default=None,
                   help="manual gain index, or 'auto' (scan, highest gain without clipping); omit = hardware AGC (not advised)")
    s.add_argument("--n", type=int, default=None, help="samples per burst (default: firmware maximum)")
    s.add_argument("--count", type=int, default=100)
    s.add_argument("--duration", type=float, default=None, help="seconds; stops at count or duration")
    s.add_argument("--per-freq", type=int, default=10, help="bursts per LO before moving to the next --freq")
    s.add_argument("--temp-every", type=int, default=25, help="read TEMP? every k bursts (0 = never)")
    s.add_argument("--out", default=str(ROOT / "data"))
    s.add_argument("--label", default="")
    s.add_argument("--note", default="", help="free text stored in session.json (attenuation, cable, position)")
    s.set_defaults(func=cmd_capture)

    s = sub.add_parser("ssb-stream", help="on-chip SSB detection (SSBDET): store only the raw I/Q around each SSB")
    common(s)
    rx(s)
    s.add_argument("--ssb", type=float, required=True, help="SSB centre in MHz, e.g. 2125.45")
    s.add_argument("--nid2", type=int, required=True, help="PSS root of the cell (PCI mod 3)")
    s.add_argument("--cfo", type=int, default=0, help="expected CFO in Hz (crystal), tracked by the firmware")
    s.add_argument("--gain", default="auto", help="manual gain index or 'auto'")
    s.add_argument("--thr", type=float, default=14.0, help="PSS threshold (noise: mean 1)")
    s.add_argument("--pre", type=int, default=20, help="us kept before the SSB")
    s.add_argument("--post", type=int, default=20, help="us kept after the SSB")
    s.add_argument("--period", type=float, default=20.0, help="SSB period of the cell in ms")
    s.add_argument("--every", type=int, default=1, help="raw cut of every k-th SSB (others header only)")
    s.add_argument("--no-track", action="store_true", help="random bursts only, no capture on the SSB grid")
    s.add_argument("--count", type=int, default=0, help="stop after this many SSB frames (0 = by duration)")
    s.add_argument("--duration", type=float, default=30.0, help="seconds")
    s.add_argument("--out", default=str(ROOT / "data"))
    s.add_argument("--label", default="")
    s.add_argument("--note", default="")
    s.set_defaults(func=cmd_ssb_stream, bits=8)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
