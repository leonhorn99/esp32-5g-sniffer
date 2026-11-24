"""Averaged spectrum of capture sessions on one absolute frequency axis (quick look).

  .venv_local/bin/python -m capture.session_spectrum data/<s1> data/<s2> ... --out results/x.png \
        [--ssb results/<s1>_ssb.csv ...]

Per session: Welch average over the bursts (1024-point Hann segments), the burst mean (DC
offset) removed first. Levels are relative to each session's median bin, because the
sessions use different gain indices; compare shapes, not levels across sessions.
SSB detections (CSV from matlab/detect_session.m) are marked at their SSB centre with the PCI.
Operator blocks are the RTR n1 assignment (2021-2044).
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from .session import load_session

N1_BLOCKS = [("A1", 2110, 2135), ("Magenta", 2135, 2150), ("Drei", 2150, 2170)]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # validated categorical slots 1-4


def session_psd(folder, nfft=1024, max_bursts=200):
    meta, recs, iq = load_session(folder)
    f_lo = recs[0]["freq_mhz"]
    fs = recs[0]["fs"]
    win = np.hanning(nfft)
    acc = np.zeros(nfft)
    for x in iq[:max_bursts]:
        x = x.astype(np.complex128)
        x = x - x.mean()
        segs = x[: (len(x) // nfft) * nfft].reshape(-1, nfft) * win
        acc += np.mean(np.abs(np.fft.fftshift(np.fft.fft(segs, axis=1), axes=1)) ** 2, axis=0)
    f = f_lo + (np.arange(nfft) - nfft // 2) * fs / nfft / 1e6
    p = 10 * np.log10(acc / np.median(acc))
    edge = np.abs(f - f_lo) > 0.9 * fs / 2e6            # drop the band edges (filter roll-off)
    return f[~edge], p[~edge], f_lo, meta["settings"].get("gain")


def read_ssb(paths):
    rows = []
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh):
                if r["sssOk"].strip().lower() in ("1", "true"):
                    rows.append((float(r["fSSB_MHz"]), int(float(r["pci"])), r["mibOk"].strip().lower() in ("1", "true")))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sessions", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ssb", nargs="*", default=[])
    ap.add_argument("--title", default="ESP32-C5 over the air, n1 downlink (measured)")
    a = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 4.2), dpi=150)
    fig.patch.set_facecolor("#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    fmin, fmax = np.inf, -np.inf
    for k, s in enumerate(a.sessions):
        f, p, f_lo, gain = session_psd(s)
        fmin, fmax = min(fmin, f.min()), max(fmax, f.max())
        ax.plot(f, p, color=COLORS[k % len(COLORS)], lw=1.2, label=f"LO {f_lo:g} MHz (gain index {gain})")
    ymin, ymax = ax.get_ylim()
    for i, (name, lo, hi) in enumerate(N1_BLOCKS):
        if hi < fmin or lo > fmax:
            continue
        ax.axvspan(max(lo, fmin), min(hi, fmax), color="#000000", alpha=0.035 if i % 2 else 0.07, lw=0)
        ax.text((max(lo, fmin) + min(hi, fmax)) / 2, ymax, name, ha="center", va="bottom", color="#555555", fontsize=9)
    seen = {}   # SSB centre -> {pci: MIB decoded at least once}
    for f_ssb, pci, mib in read_ssb(a.ssb):
        cells = seen.setdefault(round(f_ssb, 2), {})
        cells[pci] = cells.get(pci, False) or mib
    for f_ssb, cells in sorted(seen.items()):
        label = ", ".join(f"{p}{' (MIB)' if m else ''}" for p, m in sorted(cells.items()))
        ax.axvline(f_ssb, color="#333333", lw=1, ls="--")
        ax.text(f_ssb, ymin, f" SSB {f_ssb:.2f} MHz, PCI {label}", rotation=90, va="bottom", ha="right",
                color="#333333", fontsize=8)
    ax.set_xlim(fmin, fmax)
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("PSD rel. to session median (dB)")
    ax.grid(axis="y", color="#e0e0e0")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    ax.set_title(a.title, fontsize=10, pad=16)
    fig.tight_layout()
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, facecolor=fig.get_facecolor())
    print(a.out)


if __name__ == "__main__":
    main()
