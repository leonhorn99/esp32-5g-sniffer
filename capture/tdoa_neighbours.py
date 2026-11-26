"""Arrival time of neighbour cells relative to the serving cell, from the same cuts (no clock model).

  .venv_local/bin/python -m capture.tdoa_neighbours data/<SSBDET session of the serving cell> --ssb 2667.75 \
        --serving 317 --pcis 46 69 90 100 163 174 241 266 306 333 341 350 357 456 498

The neighbour list comes from the serving cell's SIB3 (capture/nr_sib1.py). In a frame-synchronised
network, a neighbour on the same frequency sends its SSB at the same time as the serving cell, so it
arrives within a few microseconds of it: inside the cut the board keeps around each serving SSB
(+-20 us at the default --pre/--post of ssb-stream, about +-15 us usable).

Per cut: the serving cell's SSB is decoded (nr_ssb.receive: refined PSS time and CFO), the cut is mixed to
SS_REF with that CFO (the board's crystal error is the same for every cell), and correlated with the SSS
of every neighbour PCI (and its PSS when its N_ID2 differs from the serving cell's: same-N_ID2 PSS are
dominated by the serving cell). Optionally the serving cell's PSS and SSS are subtracted first
(--cancel), since its cross-correlation with another SSS is only ~21 dB down. The correlation power,
normalised by the cut's power, is summed over all cuts versus the delay tau from the serving cell's
SSS; a neighbour shows up as a peak, and tau is its arrival time minus the serving cell's (times c:
the difference of the two path lengths plus any offset between the two sites' frame timing).

The statistic: peak over the median of the summed curve (outside +-1 us of the peak); about 1 for
noise. Prints one row per PCI and writes a JSON with the curves.
"""
from __future__ import annotations

import argparse
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from . import nr_ssb as nr
from .session import load_session

C = 299792458.0


@lru_cache(maxsize=256)
def templates(pci: int, fs: float) -> tuple[np.ndarray, np.ndarray]:
    """Useful parts of the PSS and SSS symbols of a PCI at fs (SS_REF at 0 Hz)."""
    w, m = nr._dft(fs, nr.SCS)
    p = (nr.pss(pci % 3) @ np.conj(w[nr.PSS_K])) / m
    s = (nr.sss(pci) @ np.conj(w[nr.PSS_K])) / m
    return p, s


def _corr(y: np.ndarray, t: np.ndarray, base: float, half: int, pw: float) -> np.ndarray | None:
    """|sum_j y[n + j] conj(t[j])|^2, normalised, for n within +-half of base; None if outside the cut."""
    i0 = int(round(base))
    lo, hi = i0 - half, i0 + half + 1
    m = len(t)
    if lo < 0 or hi - 1 + m > len(y):
        return None
    c = sliding_window_view(y[lo:hi - 1 + m], m) @ np.conj(t)
    return np.abs(c) ** 2 / (pw * np.sum(np.abs(t) ** 2))


def cut_curves(y: np.ndarray, fs: float, t_pss: float, serving: int, pcis, half: int, cancel: bool = False,
               f_ssb: float | None = None) -> dict:
    """One cut, mixed to SS_REF with the serving cell's CFO; t_pss: the serving cell's PSS useful start
    (fractional samples). Per neighbour PCI: (SSS curve, PSS curve or None) over +-half samples around
    where its SSS (PSS) lies if it arrives with the serving cell. The PSS curve only when its N_ID2
    differs (otherwise the serving cell's own PSS dominates). PCIs whose SSS window leaves the cut are
    left out."""
    if cancel:
        y = cancel_serving(y, fs, t_pss, serving, f_ssb)
    pw = float(np.mean(np.abs(y) ** 2))
    if pw <= 0:
        return {}
    ts = nr.t_sym(nr.SCS) * fs
    out = {}
    for p in pcis:
        tp, tsss = templates(int(p), float(fs))
        vs = _corr(y, tsss, t_pss + 2 * ts, half, pw)
        if vs is None:
            continue
        vp = None if p % 3 == serving % 3 else _corr(y, tp, t_pss, half, pw)
        out[int(p)] = (vs, vp)
    return out


def peak(a: np.ndarray, fs: float) -> tuple[float, float]:
    """(peak over the median outside +-1 us of the peak, delay of the peak in us from the curve's centre;
    sub-sample by a parabola through 3 samples)."""
    half = (len(a) - 1) // 2
    i = int(np.argmax(a))
    mask = np.abs(np.arange(len(a)) - i) > 1e-6 * fs
    floor = float(np.median(a[mask])) if mask.any() else float("nan")
    d = 0.0
    if 0 < i < len(a) - 1:
        den = a[i - 1] - 2 * a[i] + a[i + 1]
        d = 0.5 * (a[i - 1] - a[i + 1]) / den if den else 0.0
    return (float(a[i] / floor) if floor > 0 else float("nan")), (i + d - half) / fs * 1e6


class Accumulator:
    """Curves of many cuts summed per neighbour PCI: SSS + PSS (the delay) and SSS only (the statistic:
    specific to the PCI)."""

    def __init__(self):
        self.acc: dict[int, np.ndarray] = {}
        self.acc_s: dict[int, np.ndarray] = {}
        self.n = 0

    def add(self, curves: dict):
        if not curves:
            return
        self.n += 1
        for p, (vs, vp) in curves.items():
            if p not in self.acc_s or len(self.acc_s[p]) != len(vs):
                self.acc[p], self.acc_s[p] = np.zeros(len(vs)), np.zeros(len(vs))
            self.acc_s[p] += vs
            self.acc[p] += vs if vp is None else vs + vp

    def rows(self, fs: float) -> list[dict]:
        out = []
        for p in self.acc:
            stat, tau = peak(self.acc[p] / max(self.n, 1), fs)
            stat_s, tau_s = peak(self.acc_s[p] / max(self.n, 1), fs)
            out.append({"pci": p, "stat": stat, "tau_us": tau, "delta_m": tau * 1e-6 * C, "stat_sss": stat_s,
                        "tau_sss_us": tau_s, "curve": (self.acc[p] / max(self.n, 1)).tolist(),
                        "curve_sss": (self.acc_s[p] / max(self.n, 1)).tolist()})
        return sorted(out, key=lambda r: -r["stat_sss"])


def cancel_serving(y: np.ndarray, fs: float, t_pss: float, pci: int, f_ssb: float) -> np.ndarray:
    """y minus the serving cell's PSS and SSS (symbols 0 and 2 of its SSB), each with its
    channel measured on its own 127 subcarriers (LS, smoothed over +-4 subcarriers)."""
    w, m = nr._dft(fs, nr.SCS)
    ts = nr.t_sym(nr.SCS) * fs
    cp = nr.T_CP * fs
    out = y.copy()
    for l, seq in ((0, nr.pss(pci % 3)), (2, nr.sss(pci))):
        t_u = t_pss + l * ts                                   # useful start (fractional samples)
        n0 = int(round(t_u))
        if n0 < 0 or n0 + m > len(y):
            continue
        f = (w[nr.PSS_K] @ y[n0:n0 + m]) / m                   # the 127 subcarriers
        h = f * seq
        h = np.convolve(h, np.ones(9) / 9, mode="same")       # smoothed channel
        # regenerate the symbol including its CP and subtract it
        n = np.arange(int(round(n0 - cp)), n0 + m)
        n = n[(n >= 0) & (n < len(y))]
        k = nr.K_REL[nr.PSS_K] * nr.SCS
        sym = (h * seq) @ np.exp(2j * np.pi * np.outer(k, (n - n0)) / fs)
        out[n] -= sym
    return out


def search(folder: str | Path, f_ssb_mhz: float, serving: int, pcis: list[int], window_us: float = 15.0,
           cancel: bool = False, max_cuts: int | None = None, say=print, cut_range=None, loaded=None) -> dict:
    meta, recs, iq = loaded or load_session(folder)
    if cut_range is not None:
        recs, iq = recs[cut_range[0]:cut_range[1]], iq[cut_range[0]:cut_range[1]]
    st = meta.get("settings", {})
    f_lo = float((st.get("freq_mhz") or [0])[0]) * 1e6
    f_ssb = f_ssb_mhz * 1e6
    fs = recs[0]["fs"]
    half = int(round(window_us * 1e-6 * fs))
    acc = Accumulator()
    for r, x in list(zip(recs, iq))[:max_cuts]:
        x = np.asarray(x, complex)
        t0 = r.get("ssb_pss_in_cut")
        res = nr.receive(x, fs, t0 if t0 is not None else 0, f_ssb, f_lo, cfo=r.get("ssb_cfo_hz", 0.0),
                         nid2=serving % 3, t_first_us=0.0, decode=False)
        if res.get("pci") != serving or res.get("t_pss_us") is None:
            continue
        t_pss = res["t_pss_us"] * 1e-6 * fs
        n = np.arange(len(x))
        y = (x - x.mean()) * np.exp(-2j * np.pi * (f_ssb - f_lo + res["cfo"]) * n / fs)
        acc.add(cut_curves(y, fs, t_pss, serving, pcis, half, cancel, f_ssb))
    return {"session": str(folder), "serving": serving, "ssb_mhz": f_ssb_mhz, "cuts": acc.n, "fs": fs,
            "window_us": window_us, "cancel": cancel, "rows": acc.rows(fs)}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("session")
    p.add_argument("--ssb", type=float, required=True, help="SS_REF in MHz (serving cell and neighbours)")
    p.add_argument("--serving", type=int, required=True)
    p.add_argument("--pcis", type=int, nargs="+", required=True, help="neighbour PCIs (SIB3)")
    p.add_argument("--window", type=float, default=15.0, help="+- us around the serving cell")
    p.add_argument("--cancel", action="store_true", help="subtract the serving cell's PSS and SSS first")
    p.add_argument("--max-cuts", type=int, default=None)
    p.add_argument("--blocks", type=int, default=0, help="also split the cuts into this many blocks in time and "
                                                         "give tau per block for the detected cells (precision)")
    p.add_argument("--detect", type=float, default=8.0, help="SSS-only peak/floor that counts as detected")
    p.add_argument("--out", default=None, help="JSON (default results/tdoa_<session>.json)")
    a = p.parse_args(argv)
    r = search(a.session, a.ssb, a.serving, a.pcis, a.window, a.cancel, a.max_cuts)
    print(f"{Path(a.session).name}: {r['cuts']} cuts with PCI {a.serving}; neighbours by the SSS-only peak/floor "
          f"(about 1-3 = noise), and SSS + PSS:")
    for row in r["rows"]:
        print(f"  PCI {row['pci']:4d}: SSS {row['stat_sss']:6.1f} at {row['tau_sss_us']:+7.3f} us | SSS+PSS "
              f"{row['stat']:6.1f} at {row['tau_us']:+7.3f} us ({row['delta_m']:+8.1f} m)")
    found = [row["pci"] for row in r["rows"] if row["stat_sss"] >= a.detect]
    if a.blocks and found:
        loaded = load_session(a.session)
        n = len(loaded[1])
        edges = np.linspace(0, n, a.blocks + 1).astype(int)
        r["blocks"] = []
        print(f"per block ({n // a.blocks} cuts, about {n / a.blocks / 25:.0f} s each at 25 SSBs/s), tau in us "
              f"(SSS + PSS):")
        for k in range(a.blocks):
            b = search(a.session, a.ssb, a.serving, found, a.window, a.cancel, cut_range=(edges[k], edges[k + 1]),
                       loaded=loaded)
            taus = {row["pci"]: (row["tau_us"], row["stat"]) for row in b["rows"]}
            r["blocks"].append({"cuts": [int(edges[k]), int(edges[k + 1])], "tau_us": {p: taus[p][0] for p in found},
                                "stat": {p: taus[p][1] for p in found}})
            print(f"  block {k + 1:2d}: " + "  ".join(f"PCI {p} {taus[p][0]:+.3f} ({taus[p][1]:.0f})" for p in found))
        for p in found:
            t = np.array([blk["tau_us"][p] for blk in r["blocks"]])
            print(f"  PCI {p}: tau {t.mean():+.3f} us, std over blocks {1e3 * t.std(ddof=1):.0f} ns "
                  f"({t.std(ddof=1) * 1e-6 * C:.1f} m)")
    out = Path(a.out) if a.out else Path(__file__).resolve().parents[1] / "results" / f"tdoa_{Path(a.session).name}.json"
    out.write_text(json.dumps(r))
    print("->", out)


if __name__ == "__main__":
    main()
