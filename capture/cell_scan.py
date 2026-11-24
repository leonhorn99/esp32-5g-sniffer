"""Find NR cells: capture bursts at one or more LOs and search the whole SSB raster on the PC.

  .venv_local/bin/python -m capture.cell_scan --port /dev/ttyACM0 --freq 2131 2670
  .venv_local/bin/python -m capture.cell_scan --port /dev/ttyACM0 --freq 2131 --bursts 100 --cfo 12000 --record
  .venv_local/bin/python -m capture.cell_scan --session data/<plain burst session>      # offline

Per LO: full bursts (default 60 x 1.64 ms at 20 MS/s, about 5 SSB hits of a 20 ms cell), then
for every burst a PSS search over all GSCN positions within +-6.1 MHz of the LO (the analog
filter passes about +-8 MHz, the SSB is 3.6 MHz wide), all three N_ID2 and the CFO
hypotheses --cfo +- --cfo-span in 5 kHz steps (nr_ssb.pss_search_raster), then SSS, PBCH and
MIB (nr_ssb.receive). Prints a table of cells and writes it as JSON. 15 kHz SSBs (case A) by
default, --scs 30 --lmax 8 for 30 kHz case C (TDD bands).

The CFO is the board's crystal error times the frequency; the C5's crystal moves by several
ppm (docs/hardware.md), so the default span is +-30 kHz around --cfo.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HALF_SPAN_HZ = 6.1e6      # +-8 MHz analog passband minus half an SSB (1.8 MHz) minus margin


def gscn_number(f: float) -> int:
    if f < 3000e6:
        n = int(f // 1.2e6)
        return int(3 * n + (round((f - n * 1.2e6) / 50e3) - 3) / 2)
    return int(7499 + round((f - 3000e6) / 1.44e6))


def scan_burst(x, fs: float, f_lo: float, cfo_center: float, cfo_span: float, t0_us=None, scale: float = 1.0,
               nid2s=(0, 1, 2), scs: float = 15e3, lmax: int = 4) -> list[dict]:
    """Worker: every SSB with a PCI in one burst (15 kHz case A, or 30 kHz case C)."""
    from capture import nr_ssb as nr
    half = min(HALF_SPAN_HZ, fs / 2 - 1.9e6 * scs / 15e3)
    offs = [f - f_lo for f in nr.gscn_raster(f_lo - half, f_lo + half)]
    c0 = round(cfo_center / 250) * 250
    cfos = c0 + np.arange(-int(cfo_span // 5000), int(cfo_span // 5000) + 1) * 5000
    out = []
    for p in nr.pss_search_raster(x, fs, offs, cfos, nid2s=nid2s, scs=scs):
        f_ssb = f_lo + p["f_off"]
        r = nr.receive(x, fs, p["t_pss"], f_ssb, f_lo, cfo=p["cfo"], nid2=p["nid2"], t_first_us=t0_us, scale=scale,
                       scs=scs, lmax=lmax)
        if r.get("pci") is None:
            continue
        e = {"f_ssb": f_ssb, "pci": r["pci"], "ok": bool(r["ok"]), "snr_db": r.get("snr_db"),
             "rsrp_dbfs": r.get("rsrp_dbfs"), "cfo": r.get("cfo"), "issb": r.get("issb"), "pss_stat": p["stat"],
             "sss_stat": r.get("sss_stat"), "t_pss_us": r.get("t_pss_us")}
        if r["ok"]:
            m = r["mib"]
            e.update(sfn=m["sfn"], mib={k: m[k] for k in ("scs_common_khz", "kssb", "dmrs_typeA_pos",
                                                          "pdcch_config_sib1", "cell_barred", "intra_freq_reselection")})
        out.append(e)
    return out


def aggregate(entries: list[dict]) -> list[dict]:
    """One row per (SSB frequency, PCI), the cells with the most MIBs first."""
    cells: dict[tuple, dict] = {}
    for e in entries:
        k = (round(e["f_ssb"]), e["pci"])
        c = cells.setdefault(k, {"f_ssb_mhz": e["f_ssb"] / 1e6, "gscn": gscn_number(e["f_ssb"]), "pci": e["pci"],
                                 "ssbs": 0, "mibs": 0, "snr": [], "rsrp": [], "cfo": [], "issb": set(), "mib": None,
                                 "sfn": None})
        c["ssbs"] += 1
        c["snr"].append(e["snr_db"]); c["rsrp"].append(e["rsrp_dbfs"]); c["issb"].add(e["issb"])
        if e["ok"]:
            c["mibs"] += 1
            c["cfo"].append(e["cfo"]); c["mib"] = e["mib"]; c["sfn"] = e["sfn"]
    out = []
    for c in cells.values():
        out.append({**{k: v for k, v in c.items() if k not in ("snr", "rsrp", "cfo", "issb")},
                    "snr_db": float(np.median(c["snr"])), "rsrp_dbfs": float(np.median(c["rsrp"])),
                    "cfo_hz": float(np.median(c["cfo"])) if c["cfo"] else None, "ssb_index": sorted(c["issb"])})
    return sorted(out, key=lambda c: (-c["mibs"], -c["ssbs"]))


def scan_lo(sdr, f_lo_mhz: int, rate_code: int, n_bursts: int, pool, cfo_center: float = 0.0, cfo_span: float = 30e3,
            bits: int = 8, say=print, stop=None, writer=None, gain=None, scs: float = 15e3, lmax: int = 4) -> list[dict]:
    """Capture n_bursts full bursts at one LO (frequency and gain already set) and scan them."""
    from .ssb_live import iq_scale
    futs, entries = [], []
    for k in range(n_bursts):
        if stop and stop():
            break
        b = sdr.capture(sdr.max_samples, rate_code, bits=bits)
        if not b.crc_ok:
            continue
        if writer:
            writer.add(b, freq_mhz=f_lo_mhz, gain=gain, bandwidth_mhz=None)
        futs.append(pool.submit(scan_burst, b.iq(), b.fs, f_lo_mhz * 1e6, cfo_center, cfo_span, b.t0_us, iq_scale(b.fmt),
                                (0, 1, 2), scs, lmax))
        if k % 20 == 19:
            say(f"  LO {f_lo_mhz} MHz: {k + 1} bursts captured")
    for f in futs:
        entries += f.result()
    return aggregate(entries)


def print_table(cells, title=""):
    print(title)
    print("  SSB (MHz)   GSCN   PCI  SSBs  MIBs  SNR dB  RSRP dBFS  CFO kHz  SSB idx  MIB")
    for c in cells:
        m = c["mib"]
        mib = (f"kSSB {m['kssb']}, pdcch {m['pdcch_config_sib1']}, " + ("barred" if m["cell_barred"] else "not barred")) if m else "-"
        cfo = f"{c['cfo_hz'] / 1e3:+7.2f}" if c["cfo_hz"] is not None else "      -"
        print(f"  {c['f_ssb_mhz']:9.3f}  {c['gscn']:5d}  {c['pci']:4d}  {c['ssbs']:4d}  {c['mibs']:4d}  {c['snr_db']:+6.1f}"
              f"  {c['rsrp_dbfs']:9.1f}  {cfo}  {str(c['ssb_index']):7s}  {mib}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--port")
    p.add_argument("--freq", type=int, nargs="+", help="LOs in whole MHz")
    p.add_argument("--session", help="scan a recorded session of full bursts instead of a board")
    p.add_argument("--bursts", type=int, default=60, help="bursts per LO")
    p.add_argument("--rate", type=float, default=20, help="MS/s")
    p.add_argument("--bits", type=int, choices=(8, 10), default=8)
    p.add_argument("--gain", default="auto")
    p.add_argument("--cfo", type=float, default=0.0, help="CFO centre (Hz), e.g. the crystal error x frequency")
    p.add_argument("--ppm", type=float, default=None, help="board crystal error; sets the CFO centre per LO (-ppm x f)")
    p.add_argument("--cfo-span", type=float, default=30e3)
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--scs", type=float, default=15, choices=(15, 30), help="SSB subcarrier spacing in kHz (30: case C)")
    p.add_argument("--lmax", type=int, default=4, choices=(4, 8), help="SSBs per half frame (8: TDD above 2.4 GHz)")
    p.add_argument("--record", action="store_true", help="store the bursts as a session in data/")
    p.add_argument("--out", default=None, help="JSON with the cells (default results/cells_<time>.json)")
    a = p.parse_args(argv)
    if not a.session and not (a.port and a.freq):
        p.error("give --port and --freq, or --session")
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    from .ssb_live import _worker_init, iq_scale
    pool = ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context("spawn"), initializer=_worker_init)
    result = {"created": time.strftime("%Y-%m-%dT%H:%M:%S"), "los": []}
    t0 = time.monotonic()
    try:
        if a.session:
            from .session import load_session
            meta, recs, iq = load_session(a.session)
            by_lo: dict[float, list] = {}
            for r, x in zip(recs, iq):
                lo = float(r["freq_mhz"])
                c = a.cfo if a.ppm is None else -a.ppm * 1e-6 * lo * 1e6
                by_lo.setdefault(lo, []).append(pool.submit(scan_burst, x, r["fs"], lo * 1e6, c, a.cfo_span, r.get("t0_us"),
                                                            iq_scale(r["fmt"]), (0, 1, 2), a.scs * 1e3, a.lmax))
            for lo, futs in sorted(by_lo.items()):
                cells = aggregate([e for f in futs for e in f.result()])
                print_table(cells, f"LO {lo:g} MHz, {len(futs)} bursts:")
                result["los"].append({"lo_mhz": lo, "bursts": len(futs), "cells": cells})
        else:
            from .capture_bursts import _setup, pick_gain
            from .esp_sdr import EspSdr, rate_code_for
            from .session import SessionWriter
            sdr = EspSdr(a.port, timeout=3.0)
            device = sdr.connect()
            rate = rate_code_for(a.rate * 1e6)
            a.bandwidth, a.dccal = None, False
            writer = None
            if a.record:
                writer = SessionWriter(ROOT / "data" / f"{time.strftime('%Y%m%d_%H%M%S')}_cellscan", device,
                                       {"tool": "cell_scan", "freq_mhz": a.freq, "rate_msps": a.rate, "bits": a.bits,
                                        "gain": a.gain})
            try:
                for lo in a.freq:
                    gain = a.gain
                    if gain == "auto":
                        gain, _ = pick_gain(sdr, lo, rate, 10, a)
                    _setup(sdr, lo, a, int(gain))
                    c = a.cfo if a.ppm is None else -a.ppm * 1e-6 * lo * 1e6
                    cells = scan_lo(sdr, lo, rate, a.bursts, pool, c, a.cfo_span, a.bits, writer=writer, gain=int(gain),
                                    scs=a.scs * 1e3, lmax=a.lmax)
                    print_table(cells, f"LO {lo} MHz (gain {gain}), {a.bursts} bursts:")
                    result["los"].append({"lo_mhz": lo, "gain": int(gain), "bursts": a.bursts, "cells": cells})
            finally:
                if writer:
                    writer.close()
                sdr.close()
    finally:
        pool.shutdown()
    out = Path(a.out) if a.out else ROOT / "results" / f"cells_{time.strftime('%Y%m%d_%H%M%S')}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1))
    print(f"{sum(len(x['cells']) for x in result['los'])} cells in {time.monotonic() - t0:.0f} s -> {out}")


if __name__ == "__main__":
    main()
