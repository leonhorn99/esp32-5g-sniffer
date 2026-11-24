"""Cut the small example sessions in samples/ out of full captures (run once, from the project root).

  .venv_local/bin/python scripts/make_samples.py

The full captures stay untouched in data/ (not published). Each sample is a normal session
(session.json, bursts.bin, bursts.jsonl, and ssb_events.jsonl for SSBDET) with the original
burst numbers, plus expected_ssb.csv: the rows of the MATLAB detector's result for these
bursts (matlab/detect_session.m with the symbol-phase fix of 23 Nov 2026). The host name in
session.json is replaced.
"""
import csv
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SAMPLES = [
    # (sample folder, source session, which bursts, description)
    ("ssbdet_pci165_150cuts", "20261123_140216_ssbdet_v2", "first:150",
     "First 150 SSB cuts of an SSBDET run (on-chip detection, 20 MS/s, 8 bit) on the NR cell PCI 165, "
     "band n1, SSB at 2125.45 MHz, LO 2131 MHz, indoors with the board antenna, 23 Nov 2026."),
    ("n1_2131_20M_30bursts", "20261123_112009_n1_2131", "ssb18+12",
     "30 full bursts (32,704 samples, 20 MS/s, 10 bit) at LO 2131 MHz: the 18 that hold an SSB of PCI 165 "
     "(SSB at 2125.45 MHz) and the first 12 without one, 23 Nov 2026 (first light)."),
]


def pick(kind, recs, expected_bursts):
    if kind.startswith("first:"):
        return recs[: int(kind.split(":")[1])]
    with_ssb = [r for r in recs if r["idx"] in expected_bursts]
    without = [r for r in recs if r["idx"] not in expected_bursts][:12]
    return sorted(with_ssb + without, key=lambda r: r["idx"])


def main():
    for name, src, kind, desc in SAMPLES:
        sdir, out = ROOT / "data" / src, ROOT / "samples" / name
        with open(ROOT / "results" / f"{src}_ssb.csv") as fh:
            exp_rows = list(csv.DictReader(fh))
        exp_bursts = {int(r["burst"]) for r in exp_rows if r["sssOk"] == "1" and r["pci"] == "165"}
        recs = [json.loads(line) for line in (sdir / "bursts.jsonl").read_text().splitlines() if line.strip()]
        recs = [r for r in recs if r["crc_ok"]]
        chosen = pick(kind, recs, exp_bursts)
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True)
        meta = json.loads((sdir / "session.json").read_text())
        meta["host"] = "anonymised"
        meta["sample"] = {"source_session": src, "bursts": len(chosen), "description": desc}
        (out / "session.json").write_text(json.dumps(meta, indent=2))
        offset = 0
        with open(sdir / "bursts.bin", "rb") as fin, open(out / "bursts.bin", "wb") as fbin, \
                open(out / "bursts.jsonl", "w") as fj:
            for r in chosen:
                fin.seek(r["offset"])
                payload = fin.read(r["nbytes"])
                fbin.write(payload)
                r = dict(r, offset=offset)
                fj.write(json.dumps(r) + "\n")
                offset += len(payload)
        idx = {r["idx"] for r in chosen}
        if (sdir / "ssb_events.jsonl").exists():
            with open(out / "ssb_events.jsonl", "w") as fe:
                for line in (sdir / "ssb_events.jsonl").read_text().splitlines():
                    ev = json.loads(line)
                    if ev.get("burst_idx") in idx:
                        fe.write(line + "\n")
        with open(out / "expected_ssb.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(exp_rows[0].keys()))
            w.writeheader()
            w.writerows(r for r in exp_rows if int(r["burst"]) in idx)
        size = sum(p.stat().st_size for p in out.iterdir())
        print(f"{name}: {len(chosen)} bursts from {src}, {size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
