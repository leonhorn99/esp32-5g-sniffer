"""Stitched spectrum of band-survey sessions (capture_bursts capture with many --freq).

  .venv_local/bin/python -m capture.band_survey data/<survey1> data/<survey2> --out results/band_survey.png \
        [--ssb results/<x>_ssb.csv ...] [--threshold 6]

Per LO: Welch average over its bursts (1024-point Hann segments, burst mean removed), only the
central +-8 MHz (the analog filter rolls off beyond, and the SSB search uses no more). All LOs
of a survey use the same gain index, so levels are comparable within a survey, but not in
absolute terms (no calibration) and only roughly across frequency (the front end's response
varies). Levels are in dB over the survey's noise floor: the median, over LOs, of each LO's
10th-percentile bin.

Prints the occupied stretches (bins more than --threshold dB above the floor, merged, at least
1 MHz wide): these are where an SSB search makes sense. Band allocations in the plot are
approximate (Austria / CEPT), for orientation only.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .esp_sdr import unpack_iq

HALF_MHZ = 8.0
BANDS = [  # (label, f_lo MHz, f_hi MHz), approximate
    ("n1 DL", 2110, 2170), ("MSS (satellite) DL", 2170, 2200), ("n40 / B40 TDD", 2300, 2400),
    ("Wi-Fi 2.4 GHz", 2400, 2483.5), ("n7 / B7 UL", 2500, 2570), ("n38 / B38 TDD", 2570, 2620),
    ("n7 / B7 DL", 2620, 2690), ("radar", 2700, 2900), ("n79 (not used in Europe)", 4400, 5000),
    ("Wi-Fi 5 GHz", 5150, 5350), ("Wi-Fi 5 GHz", 5470, 5725), ("5.8 GHz SRD / Wi-Fi", 5725, 5855),
    ("ITS / C-V2X", 5855, 5925),
]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def survey_psd(folder: Path, nfft: int = 1024, max_per_lo: int = 60):
    """{LO MHz: (f_abs MHz array, mean power array)} for one survey session."""
    recs = [json.loads(line) for line in (folder / "bursts.jsonl").read_text().splitlines() if line.strip()]
    recs = [r for r in recs if r["crc_ok"]]
    win = np.hanning(nfft)
    acc: dict[float, list] = {}
    with open(folder / "bursts.bin", "rb") as fh:
        for r in recs:
            lo = float(r["freq_mhz"])
            if len(acc.get(lo, [])) >= max_per_lo:
                continue
            fh.seek(r["offset"])
            x = unpack_iq(fh.read(r["nbytes"]), r["n"], r["fmt"]).astype(complex)
            x -= x.mean()
            seg = x[: len(x) // nfft * nfft].reshape(-1, nfft) * win
            acc.setdefault(lo, []).append(np.mean(np.abs(np.fft.fftshift(np.fft.fft(seg, axis=1), axes=1)) ** 2, axis=0))
            fs = r["fs"]
    out = {}
    f_rel = (np.arange(nfft) - nfft // 2) * fs / nfft / 1e6
    keep = np.abs(f_rel) <= HALF_MHZ
    keep &= np.abs(f_rel) > 0.1          # DC bin (residual LO leakage)
    for lo, ps in acc.items():
        out[lo] = (lo + f_rel[keep], np.mean(ps, axis=0)[keep])
    return out


def occupied(f: np.ndarray, db: np.ndarray, thr: float, min_mhz: float = 1.0):
    order = np.argsort(f)
    f, db = f[order], db[order]
    hot = db > thr
    segs, start = [], None
    for i, h in enumerate(hot):
        if h and start is None:
            start = i
        if (not h or i == len(hot) - 1) and start is not None:
            end = i if not h else i + 1
            if f[end - 1] - f[start] >= min_mhz:
                segs.append((f[start], f[end - 1], float(np.max(db[start:end])), float(np.median(db[start:end]))))
            start = None
    return segs


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("sessions", nargs="+")
    p.add_argument("--out", default="results/band_survey.png")
    p.add_argument("--ssb", nargs="*", default=[], help="detection CSVs (matlab/detect_session.m) to mark")
    p.add_argument("--threshold", type=float, default=6.0, help="dB over the floor for 'occupied'")
    a = p.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surveys = []
    for s in a.sessions:
        psd = survey_psd(Path(s))
        floor = np.median([np.percentile(v[1], 10) for v in psd.values()])
        f = np.concatenate([v[0] for v in psd.values()])
        db = np.concatenate([10 * np.log10(v[1] / floor) for v in psd.values()])
        surveys.append((Path(s).name, psd, floor, f, db))
        print(f"{Path(s).name}: {len(psd)} LOs {min(psd):.0f}-{max(psd):.0f} MHz")
        for f0, f1, mx, md in occupied(f, db, a.threshold):
            print(f"  occupied {f0:8.2f} - {f1:8.2f} MHz ({f1 - f0:5.1f} MHz), median {md:4.1f} dB, max {mx:4.1f} dB over the floor")

    ssb = []
    for path in a.ssb:
        with open(path) as fh:
            for r in csv.DictReader(fh):
                if r.get("sssOk") in ("1", "true") and r.get("pci") not in ("", "NaN"):
                    ssb.append((float(r["fSSB_MHz"]), int(float(r["pci"]))))
    ssb = sorted(set(ssb))

    fig, axs = plt.subplots(len(surveys), 1, figsize=(11, 3.0 * len(surveys)), squeeze=False)
    for ax, (name, psd, floor, f, db) in zip(axs[:, 0], surveys):
        lo_all = sorted(psd)
        for k, lo in enumerate(lo_all):
            fr, pw = psd[lo]
            order = np.argsort(fr)
            ax.plot(fr[order], 10 * np.log10(pw[order] / floor), lw=0.7, color=COLORS[k % 2])
        x0, x1 = min(lo_all) - HALF_MHZ, max(lo_all) + HALF_MHZ
        ymax = max(20, float(np.percentile(db, 99.9)) + 3)
        for k, (lab, b0, b1) in enumerate(BANDS):
            if b1 < x0 or b0 > x1:
                continue
            ax.axvspan(max(b0, x0), min(b1, x1), color="#eef3fb" if k % 2 else "#f6f1e8", zorder=0)
            ax.text((max(b0, x0) + min(b1, x1)) / 2, ymax * 0.97, lab, ha="center", va="top", fontsize=7, color="#555")
        for fs_mhz, pci in ssb:
            if x0 <= fs_mhz <= x1:
                ax.axvline(fs_mhz, color="#c00", lw=0.8, ls="--")
                ax.text(fs_mhz, ymax * 0.8, f"PCI {pci}", rotation=90, fontsize=7, color="#c00", ha="right")
        ax.axhline(a.threshold, color="#999", lw=0.6, ls=":")
        ax.set_xlim(x0, x1)
        ax.set_ylim(-3, ymax)
        ax.set_ylabel("dB over noise floor")
        ax.set_title(f"{name}: {len(lo_all)} LOs, central ±{HALF_MHZ:.0f} MHz each (colours alternate per LO)",
                     fontsize=9, loc="left")
        ax.grid(True, color="#eee", lw=0.6)
    axs[-1, 0].set_xlabel("frequency (MHz)")
    fig.tight_layout()
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=130)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
