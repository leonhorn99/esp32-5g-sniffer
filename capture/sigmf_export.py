"""Export a capture session as a SigMF recording (https://sigmf.org), for other SDR tools.

  .venv_local/bin/python -m capture.sigmf_export samples/ssbdet_pci165_150cuts results/ssb_sample
  .venv_local/bin/python -m capture.sigmf_export data/<session> results/<name> --ssb results/<session>_ssb.csv

Writes <name>.sigmf-data and <name>.sigmf-meta:
- the data file holds all bursts or cuts back to back: CAP16 as ci8, CAP20 as ci16_le
  (10-bit values), I then Q, Q already in the true orientation;
- one SigMF capture segment per burst, with its centre frequency (the LO, or the SSB centre
  for decimated cuts) and the board's timestamp (extension fields esp32:t0_us and
  esp32:burst; core:datetime is the host's receive time, so only approximate);
- one annotation per SSB: from the SSBDET metadata (PSS position in the cut) and, with --ssb,
  from a detection CSV (matlab/detect_session.m, or samples/*/expected_ssb.csv) with the PCI,
  SFN and SSB index as the label.
The sessions' bursts are not contiguous in time; SigMF readers that assume one continuous
stream see them as one, so use the capture segments.
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from .session import load_session

SSB_S = 4 * 2192 / 30.72e6       # four OFDM symbols at 15 kHz with normal CP (285.4 us)
CP_S = 144 / 30.72e6


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("session")
    p.add_argument("out", help="output path without extension")
    p.add_argument("--ssb", help="detection CSV for labelled SSB annotations")
    a = p.parse_args(argv)
    meta, recs, iq = load_session(a.session)
    if not recs:
        raise SystemExit("empty session")
    fs = {r["fs"] for r in recs}
    if len(fs) != 1:
        raise SystemExit(f"bursts at several sample rates {fs}: SigMF needs one per recording")
    fs = fs.pop()
    fmt = recs[0]["fmt"]
    dt, dtype = (np.int8, "ci8") if fmt == "CAP16" else (np.dtype("<i2"), "ci16_le")
    st = meta.get("settings", {})
    f_ssb = st.get("ssb_mhz")
    det = {}
    if a.ssb:
        with open(a.ssb) as fh:
            for r in csv.DictReader(fh):
                if r.get("sssOk") in ("1", "true"):
                    det.setdefault(int(r["burst"]), []).append(r)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    captures, annotations, start = [], [], 0
    with open(out.with_suffix(".sigmf-data"), "wb") as fd:
        for rec, x in zip(recs, iq):
            inter = np.empty(2 * len(x), dtype=dt)
            inter[0::2] = np.real(x).astype(dt)
            inter[1::2] = np.imag(x).astype(dt)
            fd.write(inter.tobytes())
            dec = (rec.get("ssb_R") or 1) > 1
            fc = f_ssb * 1e6 if dec and f_ssb else float(rec["freq_mhz"]) * 1e6
            cap = {"core:sample_start": start, "core:frequency": fc, "esp32:burst": rec["idx"],
                   "esp32:t0_us": rec.get("t0_us"), "esp32:gain_index": rec.get("gain")}
            if rec.get("host_time"):
                cap["core:datetime"] = datetime.fromtimestamp(rec["host_time"], timezone.utc).isoformat().replace("+00:00", "Z")
            captures.append(cap)
            ssbs = []
            if rec.get("ssb_pss_in_cut") is not None:
                ssbs.append((float(rec["ssb_pss_in_cut"]), None))
            for r in det.get(rec["idx"], []):
                if rec.get("t0_us") is not None:
                    ssbs = [s for s in ssbs if abs(s[0] - (float(r["tPss_us"]) - rec["t0_us"]) * fs / 1e6) > 5] + \
                        [((float(r["tPss_us"]) - rec["t0_us"]) * fs / 1e6, r)]
            for pos, r in ssbs:
                s0 = int(round(pos - CP_S * fs))
                if s0 < 0 or s0 >= len(x):
                    continue
                ann = {"core:sample_start": start + s0, "core:sample_count": int(round(SSB_S * fs)),
                       "core:label": "NR SSB"}
                fcen = f_ssb * 1e6 if f_ssb else None
                if r is not None:
                    fcen = float(r["fSSB_MHz"]) * 1e6
                    lab = f"NR SSB PCI {r['pci']} SSB {r['iSSB']}"
                    if r.get("mibOk") in ("1", "true"):
                        lab += f" SFN {int(float(r['sfn']))}"
                    ann["core:label"] = lab
                if fcen:
                    ann["core:freq_lower_edge"] = round(fcen - 1.8e6)
                    ann["core:freq_upper_edge"] = round(fcen + 1.8e6)
                annotations.append(ann)
            start += len(x)
    dev = meta.get("device", {})
    m = {"global": {"core:datatype": dtype, "core:sample_rate": fs, "core:version": "1.2.0",
                    "core:description": f"ESP32-C5 capture session {Path(a.session).name}"
                                        + (" (decimated SSB cuts, SSB centre at 0 Hz)" if (recs[0].get("ssb_R") or 1) > 1 else ""),
                    "core:hw": f"ESP32-C5, {dev.get('info', '')} {dev.get('fwopt', '')}".strip(),
                    "core:recorder": "esp32-5g-sniffer", "core:num_channels": 1,
                    "core:extensions": [{"name": "esp32", "version": "0.1.0", "optional": True}],
                    "esp32:settings": st},
         "captures": captures, "annotations": sorted(annotations, key=lambda q: q["core:sample_start"])}
    out.with_suffix(".sigmf-meta").write_text(json.dumps(m, indent=1))
    print(f"wrote {out}.sigmf-data/.sigmf-meta: {start} samples ({dtype}, {fs / 1e6:g} MS/s), "
          f"{len(captures)} captures, {len(annotations)} SSB annotations")


if __name__ == "__main__":
    main()
