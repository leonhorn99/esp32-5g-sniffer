"""NR SS/PBCH block receiver in NumPy (15 kHz, case A, L_max 4: band n1 and similar).

The same chain as matlab/ssb_detect.m, for live use without MATLAB (capture/ssb_live.py):
PSS -> SSS (PCI) -> PBCH DM-RS (SSB index, half frame) -> channel estimate -> PBCH ->
polar list decoding -> BCH -> MIB. Everything works on one cut of raw I/Q around a known
PSS position (an SSBDET frame), or on a whole burst after `pss_search`.

Conventions:
- SSB grid: 240 subcarriers x 4 symbols, k = 0..239, SS_REF (the SSB "centre") = k 120.
- TS 38.211 5.4 symbol phase: every subcarrier's phase restarts at the useful start of every
  symbol at its absolute RF frequency, so a receiver tuned to SS_REF sees a phase step of
  2 pi f_SSB T_sym from symbol to symbol. It is removed here with the absolute SSB frequency
  (without it, the PSS->SSS CFO refinement is off by -f_SSB mod 7.0 kHz: +3.93 kHz at
  2125.45 MHz, found on 23 Nov 2026).
- Bits and sequences follow TS 38.211 / 38.212 and are checked against the 5G Toolbox
  (tests/test_nr_ssb.py, reference vectors from matlab/export_pbch_tables.m).
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from . import nr_pbch_tables as _T

SCS = 15e3
T_U = 1 / SCS                       # useful symbol length
T_SYM = 2192 / 30.72e6              # normal-CP symbol (all four SSB symbols in case A)
T_CP = 144 / 30.72e6
BACKOFF = 9 / 7.68e6                # FFT window starts this far into the CP (as ssb_detect.m)
PSS_K = np.arange(56, 183)          # PSS / SSS subcarriers
K_REL = np.arange(240) - 120        # subcarrier offset from SS_REF
LMAX = 4
G_TAB = (16, 23, 18, 17, 8, 30, 10, 6, 24, 7, 0, 5, 3, 2, 1, 4, 9, 11, 12, 13, 14, 15, 19, 20, 21,
         22, 25, 26, 27, 28, 29, 31)  # TS 38.212 Table 7.1.1-1
CRC24C = 0xB2B117
PSS_THRESHOLD = 18.0                # as ssb_detect.m: unit-mean exponential under noise
SSS_THRESHOLD = 13.0


# ----------------------------------------------------------------------------- sequences
@lru_cache(maxsize=4096)
def gold(c_init: int, n: int) -> np.ndarray:
    """Pseudo-random sequence c(0..n-1) of TS 38.211 5.2.1."""
    nc = 1600
    x1 = np.zeros(nc + n + 31, np.uint8)
    x2 = np.zeros(nc + n + 31, np.uint8)
    x1[0] = 1
    x2[:31] = [(c_init >> i) & 1 for i in range(31)]
    for i in range(nc + n):
        x1[i + 31] = x1[i + 3] ^ x1[i]
        x2[i + 31] = x2[i + 3] ^ x2[i + 2] ^ x2[i + 1] ^ x2[i]
    out = x1[nc:nc + n] ^ x2[nc:nc + n]
    out.setflags(write=False)
    return out


def _mseq(taps: tuple[int, ...], init: tuple[int, ...]) -> np.ndarray:
    x = np.zeros(127, np.int8)
    x[:7] = init
    for i in range(120):
        x[i + 7] = sum(x[i + t] for t in taps) % 2
    return x


_X_PSS = _mseq((4, 0), (0, 1, 1, 0, 1, 1, 1))
_X0_SSS = _mseq((4, 0), (1, 0, 0, 0, 0, 0, 0))
_X1_SSS = _mseq((1, 0), (1, 0, 0, 0, 0, 0, 0))


def pss(nid2: int) -> np.ndarray:
    n = np.arange(127)
    return 1.0 - 2.0 * _X_PSS[(n + 43 * nid2) % 127]


def sss(pci: int) -> np.ndarray:
    nid1, nid2 = divmod(pci, 3)
    m0 = 15 * (nid1 // 112) + 5 * nid2
    m1 = nid1 % 112
    n = np.arange(127)
    return (1.0 - 2.0 * _X0_SSS[(n + m0) % 127]) * (1.0 - 2.0 * _X1_SSS[(n + m1) % 127])


@lru_cache(maxsize=3)
def sss_matrix(nid2: int) -> np.ndarray:
    """(336, 127): row nid1 = SSS of PCI 3 nid1 + nid2."""
    return np.array([sss(3 * nid1 + nid2) for nid1 in range(336)])


@lru_cache(maxsize=8192)
def pbch_dmrs(pci: int, ibar: int) -> np.ndarray:
    """144 DM-RS symbols (TS 38.211 7.4.1.4.1); ibar = i_SSB + 4 n_hf for L_max 4."""
    c_init = (1 << 11) * (ibar + 1) * (pci // 4 + 1) + (1 << 6) * (ibar + 1) + (pci % 4)
    c = gold(c_init, 288).astype(float)
    return ((1 - 2 * c[0::2]) + 1j * (1 - 2 * c[1::2])) / np.sqrt(2)


@lru_cache(maxsize=4)
def ssb_indices(v: int) -> tuple[np.ndarray, np.ndarray]:
    """(PBCH, DM-RS) indices into the 240 x 4 grid, 0-based, column-major (k + 240 l),
    in mapping order (k first, then l). v = PCI mod 4."""
    pbch, dmrs = [], []
    for l in (1, 2, 3):
        for k in range(240):
            if l == 2 and 48 <= k < 192:
                continue
            (dmrs if k % 4 == v else pbch).append(k + 240 * l)
    return np.array(pbch), np.array(dmrs)


def qpsk(bits: np.ndarray) -> np.ndarray:
    b = np.asarray(bits, float)
    return ((1 - 2 * b[0::2]) + 1j * (1 - 2 * b[1::2])) / np.sqrt(2)


# ----------------------------------------------------------------------------- CRC, polar
def crc24c(bits) -> np.ndarray:
    reg = 0
    for b in bits:
        fb = ((reg >> 23) & 1) ^ int(b)
        reg = (reg << 1) & 0xFFFFFF
        if fb:
            reg ^= CRC24C
    return np.array([(reg >> (23 - i)) & 1 for i in range(24)], np.uint8)


_INFO = np.array(_T.INFO_POS)
_RM = np.array(_T.RM_MAP)
_FROZEN = np.ones(_T.N, bool)
_FROZEN[_INFO] = False


def polar_transform(u: np.ndarray) -> np.ndarray:
    """x = u G_N (G_N = [[1,0],[1,1]] Kronecker power, TS 38.212 5.3.1); its own inverse."""
    x = np.array(u, np.uint8)
    lead = x.shape[:-1]
    n = x.shape[-1]
    s = 1
    while s < n:
        y = x.reshape(lead + (n // (2 * s), 2, s))
        y[..., 0, :] ^= y[..., 1, :]
        s *= 2
    return x


def _scl(alpha, pm, frozen, lmax):
    """Successive-cancellation list decoding of one subtree (min-sum, natural order).
    alpha (P, n) LLRs of the subtree's code bits for P paths; returns the re-encoded code
    bits (P', n), the path metrics (P',) and each surviving path's parent index."""
    p, n = alpha.shape
    if frozen.all():                                   # rate-0 subtree: all bits 0
        return np.zeros((p, n), np.uint8), pm + np.sum(np.where(alpha < 0, -alpha, 0.0), axis=1), np.arange(p)
    if n == 1:                                         # information bit: split every path
        a = alpha[:, 0]
        allpm = np.concatenate([pm + np.where(a < 0, -a, 0.0), pm + np.where(a > 0, a, 0.0)])
        keep = np.argsort(allpm, kind="stable")[:min(2 * p, lmax)]
        return (keep >= p).astype(np.uint8)[:, None], allpm[keep], keep % p
    h = n // 2
    a1, a2 = alpha[:, :h], alpha[:, h:]
    b1, pm, p1 = _scl(np.sign(a1) * np.sign(a2) * np.minimum(np.abs(a1), np.abs(a2)), pm, frozen[:h], lmax)
    a1, a2 = a1[p1], a2[p1]
    b2, pm, p2 = _scl(a2 + (1.0 - 2.0 * b1) * a1, pm, frozen[h:], lmax)
    return np.concatenate([b1[p2] ^ b2, b2], axis=1), pm, p1[p2]


# ----------------------------------------------------------------------------- BCH
_SFN_POS = (1, 2, 3, 4, 5, 6, 24, 25, 26, 27)       # SFN bits in the 32-bit payload, MSB first
_UNSCRAMBLED = (G_TAB[10], G_TAB[7], G_TAB[8])      # half frame, 3rd and 2nd LSB of SFN


def _interleave_map() -> np.ndarray:
    """dest[i]: position of payload bit i after the PBCH payload interleaving (7.1.1)."""
    dest = np.zeros(32, int)
    j_sfn, j_hrf, j_ssb, j_other = 0, 10, 11, 14
    for i in range(32):
        if i in _SFN_POS:
            dest[i] = G_TAB[j_sfn]
            j_sfn += 1
        elif i == 28:
            dest[i] = G_TAB[j_hrf]
        elif 29 <= i <= 31:
            dest[i] = G_TAB[j_ssb]
            j_ssb += 1
        else:
            dest[i] = G_TAB[j_other]
            j_other += 1
    return dest


_DEST = _interleave_map()


def _payload_scrambling(pci: int, v: int) -> np.ndarray:
    c = gold(pci, 29 * 4)
    s = np.zeros(32, np.uint8)
    j = 0
    for i in range(32):
        if i not in _UNSCRAMBLED:
            s[i] = c[j + 29 * v]
            j += 1
    return s


def bch_encode(mib24, sfn: int, hrf: int, kssb: int, pci: int) -> np.ndarray:
    """864 coded BCH bits (TS 38.212 7.1, L_max 4), as nrBCH(mib24, sfn, hrf, 4, kssb, pci)."""
    a = np.concatenate([np.asarray(mib24, np.uint8),
                        [(sfn >> 3) & 1, (sfn >> 2) & 1, (sfn >> 1) & 1, sfn & 1, hrf, (kssb >> 4) & 1, 0, 0]])
    ap = np.zeros(32, np.uint8)
    ap[_DEST] = a
    ap ^= _payload_scrambling(pci, (sfn >> 1) & 3)
    c = np.concatenate([ap, crc24c(ap)])
    u = np.zeros(_T.N, np.uint8)
    u[_INFO] = c
    return polar_transform(u)[_RM]


def pbch_scrambling(pci: int, v: int) -> np.ndarray:
    """TS 38.211 7.3.3.1, L_max 4: v = i_SSB mod 4."""
    return gold(pci, 864 * 4)[864 * v:864 * (v + 1)]


def bch_decode(llr: np.ndarray, pci: int, list_length: int = 8) -> dict | None:
    """Rate recovery, CRC-aided polar list decoding, descrambling and deinterleaving of
    the 864 PBCH LLRs (after PBCH descrambling; positive = bit 0). None if no path
    passes the CRC."""
    y = np.bincount(_RM, weights=llr, minlength=_T.N)
    x, pm, _ = _scl(y[None, :], np.zeros(1), _FROZEN, list_length)
    for i in np.argsort(pm, kind="stable"):
        c = polar_transform(x[i])[_INFO]
        if np.any(crc24c(c[:32]) != c[32:]):
            continue
        ap = c[:32].copy()
        v = 2 * int(ap[G_TAB[7]]) + int(ap[G_TAB[8]])
        ap ^= _payload_scrambling(pci, v)
        a = ap[_DEST]
        sfn = int("".join(str(int(b)) for b in a[list(_SFN_POS)]), 2)
        return parse_mib(a[:24], sfn, int(a[28]), int(a[29]), path_metric=float(pm[i]))
    return None


def parse_mib(mib24, sfn: int, hrf: int, kssb_msb: int, path_metric: float = float("nan")) -> dict:
    """MIB fields (TS 38.331, UPER: choice bit, then the MIB)."""
    b = [int(x) for x in mib24]
    val = lambda i, n: int("".join(map(str, b[i:i + n])), 2)  # noqa: E731
    pdcch = val(13, 8)
    return {"sfn": sfn, "hrf": hrf, "scs_common_khz": (15, 30)[b[7]],
            "kssb": val(8, 4) + 16 * kssb_msb, "dmrs_typeA_pos": (2, 3)[b[12]],
            "pdcch_config_sib1": pdcch, "coreset0": pdcch >> 4, "searchspace0": pdcch & 15,
            "cell_barred": b[21] == 0, "intra_freq_reselection": ("allowed", "notAllowed")[b[22]],
            "bits": "".join(map(str, b)), "path_metric": path_metric}


# ----------------------------------------------------------------------------- timing
def useful_start_s(symbol: int) -> float:
    """Start of the useful part of OFDM symbol `symbol` (counted from the half-frame start), 15 kHz."""
    slot, l = divmod(symbol, 14)
    cp = lambda m: 160 if m in (0, 7) else 144  # noqa: E731
    return slot * 1e-3 + (sum(2048 + cp(m) for m in range(l)) + cp(l)) / 30.72e6


SSB_FIRST_SYMBOL = (2, 8, 16, 22)   # case A, L_max 4


# ----------------------------------------------------------------------------- front end
@lru_cache(maxsize=8)
def _dft(fs: float) -> tuple[np.ndarray, int]:
    m = int(round(T_U * fs))
    w = np.exp(-2j * np.pi * np.outer(K_REL * SCS, np.arange(m)) / fs)
    return w, m


@lru_cache(maxsize=8)
def _pss_time(nid2: int, fs: float) -> np.ndarray:
    """Useful part of the PSS symbol at sample rate fs (SS_REF at 0 Hz)."""
    w, m = _dft(fs)
    return (pss(nid2) @ np.conj(w[PSS_K])) / m


def grid(y: np.ndarray, fs: float, t_pss: float, f_rf: float) -> np.ndarray:
    """240 x 4 SSB grid from y (SS_REF at 0 Hz, CFO removed); t_pss = PSS useful start in
    samples (fractional). f_rf: absolute SS_REF frequency for the 38.211 symbol phase."""
    w, m = _dft(fs)
    g = np.zeros((240, 4), complex)
    for l in range(4):
        t_us = t_pss + l * T_SYM * fs
        n0 = int(round(t_us - BACKOFF * fs))
        if n0 < 0 or n0 + m > len(y):
            raise ValueError("SSB not inside the cut")
        f = (w @ y[n0:n0 + m]) / m
        g[:, l] = f * np.exp(2j * np.pi * K_REL * SCS * (t_us - n0) / fs) \
            * np.exp(2j * np.pi * ((f_rf * l * T_SYM) % 1.0))
    return g


def _smooth(vals: np.ndarray, ks: np.ndarray, half: int) -> tuple[np.ndarray, np.ndarray]:
    """Average of pilot values within +-half subcarriers, for every k = 0..239."""
    s = np.zeros(240 + 2 * half + 1, complex)
    c = np.zeros(240 + 2 * half + 1)
    np.add.at(s, ks + half + 1, vals)
    np.add.at(c, ks + half + 1, 1)
    cs, cc = np.cumsum(s), np.cumsum(c)
    k = np.arange(240)
    tot = cs[k + 2 * half + 1] - cs[k]
    cnt = cc[k + 2 * half + 1] - cc[k]
    return np.where(cnt > 0, tot / np.maximum(cnt, 1), 0), cnt


def pss_search(x: np.ndarray, fs: float, offset_hz: float, nid2s=(0, 1, 2), cfo_max: float = 30e3,
               cfo_step: float = 5e3, threshold: float = PSS_THRESHOLD, max_peaks: int = 6) -> list[dict]:
    """PSS peaks in a whole burst: FFT correlation over N_ID2 and CFO hypotheses after a
    +-1.1 MHz brick-wall low-pass, statistic normalised to unit mean under white noise
    (as ssb_detect.m). offset_hz = SS_REF - LO. Returns dicts t_pss (samples), nid2, stat, cfo."""
    x = np.asarray(x, complex) - np.mean(x)
    n = np.arange(len(x))
    y = x * np.exp(-2j * np.pi * offset_hz * n / fs)
    w, m = _dft(fs)
    nf = 1 << int(np.ceil(np.log2(len(y) + m)))
    f = np.fft.fftfreq(nf, 1 / fs)
    lp = np.abs(f) <= 1.1e6
    beta = 2.2e6 / fs
    peaks = []
    best = {}
    for cfo in np.arange(-cfo_max, cfo_max + 1, cfo_step):
        z = y * np.exp(-2j * np.pi * cfo * n / fs)
        zf = np.fft.fft(z, nf)
        zl = np.fft.ifft(zf * lp)[:len(y)]
        e = np.convolve(np.abs(zl) ** 2, np.ones(m), "valid")
        for nid2 in nid2s:
            r = _pss_time(nid2, fs)
            c = np.fft.ifft(zf * np.conj(np.fft.fft(r, nf)))[:len(e)]
            st = np.abs(c) ** 2 * m * beta / (np.maximum(e, 1e-12) * np.sum(np.abs(r) ** 2))
            if nid2 not in best:
                best[nid2] = (st, np.full(len(st), cfo))
            else:
                b, bc = best[nid2]
                upd = st > b
                b[upd] = st[upd]
                bc[upd] = cfo
    sym = int(T_SYM * fs)
    for nid2, (st, cf) in best.items():
        st = st.copy()
        for _ in range(max_peaks):
            i = int(np.argmax(st))
            if st[i] < threshold:
                break
            peaks.append({"t_pss": float(i), "nid2": nid2, "stat": float(st[i]), "cfo": float(cf[i])})
            st[max(0, i - sym):i + sym] = 0
    return sorted(peaks, key=lambda p: p["t_pss"])


def _refine_peak(y: np.ndarray, r: np.ndarray, t: float, fs: float, half_us: float = 0.5) -> float:
    m = len(r)
    w = max(1, int(round(half_us * 1e-6 * fs)))
    lags = np.arange(int(round(t)) - w, int(round(t)) + w + 1)
    lags = lags[(lags >= 0) & (lags + m <= len(y))]
    if len(lags) < 3:
        return t
    c = np.array([abs(np.vdot(r, y[i:i + m])) ** 2 for i in lags])
    i = int(np.argmax(c))
    if 0 < i < len(c) - 1:
        d = c[i - 1] - 2 * c[i] + c[i + 1]
        return float(lags[i] + (0.5 * (c[i - 1] - c[i + 1]) / d if d < 0 else 0.0))
    return float(lags[i])


def receive(iq: np.ndarray, fs: float, t_pss: float, f_ssb: float, f_lo: float, cfo: float = 0.0,
            nid2: int | None = None, t_first_us: float | None = None, decode: bool = True,
            list_length: int = 8) -> dict:
    """Process one SSB. iq: raw samples; t_pss: PSS useful start (samples, from the
    firmware or pss_search); cfo: CFO estimate (Hz) to start from. Returns a dict with
    the grid, channel, PCI, SSB index, CFO, SNR, MIB (or None) and timing."""
    x = np.asarray(iq, complex)
    x = x - np.mean(x)
    n = np.arange(len(x))
    f_shift = f_ssb - f_lo
    out = {"ok": False, "fs": fs, "f_ssb": f_ssb, "cfo_in": cfo, "nid2": nid2, "pci": None, "mib": None}

    def mixed(c):
        return x * np.exp(-2j * np.pi * (f_shift + c) * n / fs)

    # CFO inside the PSS symbol (two halves: free of the 38.211 symbol phase)
    y = mixed(cfo)
    w, m = _dft(fs)
    p0 = int(round(t_pss))
    if p0 < 0 or p0 + m > len(y):
        out["error"] = "PSS outside the cut"
        return out
    if nid2 is None:
        nid2 = int(np.argmax([abs(np.vdot(_pss_time(k, fs), y[p0:p0 + m])) for k in range(3)]))
        out["nid2"] = nid2
    r = _pss_time(nid2, fs)
    h = m // 2
    c1 = np.vdot(r[:h], y[p0:p0 + h])
    c2 = np.vdot(r[h:], y[p0 + h:p0 + m])
    cfo += np.angle(c2 * np.conj(c1)) / (2 * np.pi * h / fs)
    y = mixed(cfo)
    # fine timing: PSS correlation peak within +-0.5 us at the raw rate, parabolic interpolation
    # (the phase slope of the PSS channel was noisier: 0.51 us RMS in the SFN-anchored frame
    # timing against 0.29 us for this, session 20261123_140216_ssbdet_v2)
    t_pss = _refine_peak(y, r, t_pss, fs)
    try:
        g = grid(y, fs, t_pss, f_ssb)
    except ValueError as e:
        out["error"] = str(e)
        return out
    h0 = g[PSS_K, 0] * pss(nid2)
    # SSS -> PCI, against the smoothed PSS channel
    hs = np.convolve(h0, np.ones(9) / 9, "same")
    a = g[PSS_K, 2] * np.conj(hs)
    cs = sss_matrix(nid2) @ a
    st = np.abs(cs) ** 2 / max(np.sum(np.abs(a) ** 2), 1e-30)
    order = np.argsort(st)[::-1]
    out.update(sss_stat=float(st[order[0]]), sss_ratio=float(st[order[0]] / st[order[1]]))
    out["pss_power"] = float(np.mean(np.abs(h0) ** 2))
    if st[order[0]] < SSS_THRESHOLD:
        out.update(grid=g, cfo=cfo, t_pss=t_pss, error="SSS below threshold")
        return out
    pci = 3 * int(order[0]) + nid2
    out["pci"] = pci
    # CFO refinement from the PSS -> SSS phase (two symbols; unambiguous within +-3.5 kHz)
    cfo += np.angle(cs[order[0]]) / (2 * np.pi * 2 * T_SYM)
    y = mixed(cfo)
    g = grid(y, fs, t_pss, f_ssb)
    out.update(grid=g, cfo=float(cfo), t_pss=float(t_pss))
    gf = g.flatten(order="F")
    pbch_i, dmrs_i = ssb_indices(pci % 4)
    dk, dl = dmrs_i % 240, dmrs_i // 240
    # SSB index and half frame from the DM-RS: coherent within 12-subcarrier chunks over
    # the three symbols (the channel is static over the SSB once CFO and symbol phase are gone)
    metric = np.zeros(8)
    chunk = dk // 12
    for ibar in range(8):
        v = gf[dmrs_i] * np.conj(pbch_dmrs(pci, ibar))
        s = np.zeros(20, complex)
        np.add.at(s, chunk, v)
        metric[ibar] = np.sum(np.abs(s) ** 2) / np.sum(np.abs(v) ** 2)
    cand = np.argsort(metric)[::-1]
    out["dmrs_metric"] = metric
    # channel estimate from all pilots (PSS, SSS, DM-RS of the best candidate)
    best = None
    for ibar in cand[:3] if decode else cand[:1]:
        ref = pbch_dmrs(pci, int(ibar))
        pil_k = np.concatenate([PSS_K, PSS_K, dk])
        pil_v = np.concatenate([g[PSS_K, 0] * pss(nid2), g[PSS_K, 2] * sss(pci), gf[dmrs_i] * np.conj(ref)])
        hk, cnt = _smooth(pil_v, pil_k, 6)
        resid = gf[dmrs_i] - hk[dk] * ref
        n_eff = np.mean(cnt[dk])
        nvar = max(np.mean(np.abs(resid) ** 2) * n_eff / max(n_eff - 1, 1), 1e-30)
        pk, pl = pbch_i % 240, pbch_i // 240
        hp = hk[pk]
        yp = gf[pbch_i]
        snr = float(np.mean(np.abs(hp) ** 2) / nvar)
        z = yp * np.conj(hp)
        llr = np.empty(864)
        llr[0::2] = 2 * np.sqrt(2) * z.real / nvar
        llr[1::2] = 2 * np.sqrt(2) * z.imag / nvar
        # equalised symbols for display: MMSE (weak subcarriers shrink towards 0 instead of
        # being blown up by 1/H), scaled by one factor so the average symbol is unbiased
        w = np.abs(hp) ** 2 / (np.abs(hp) ** 2 + nvar)
        eq = z / (np.abs(hp) ** 2 + nvar) / max(np.mean(w), 1e-9)
        res = {"ibar": int(ibar), "snr_db": 10 * np.log10(snr), "channel": hk, "eq": eq, "nvar": nvar}
        if best is None:
            best = res
        if not decode:
            break
        v = int(ibar) % 4
        mib = bch_decode(llr * (1 - 2.0 * pbch_scrambling(pci, v)), pci, list_length)
        if mib is not None:
            best = res
            best["mib"] = mib
            break
    out.update({k: best[k] for k in ("ibar", "snr_db", "channel", "eq", "nvar")})
    out["issb"] = out["ibar"] % 4
    out["hrf_dmrs"] = out["ibar"] // 4
    out["mib"] = best.get("mib")
    out["ok"] = out["mib"] is not None
    hrf = out["mib"]["hrf"] if out["ok"] else out["hrf_dmrs"]
    if t_first_us is not None:
        t_pss_us = t_first_us + t_pss / fs * 1e6
        out["t_pss_us"] = t_pss_us
        out["t_frame_us"] = t_pss_us - 1e6 * useful_start_s(SSB_FIRST_SYMBOL[out["issb"]]) - 5000 * hrf
    return out
