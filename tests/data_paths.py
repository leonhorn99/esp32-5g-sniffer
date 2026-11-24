"""Where the tests find recorded captures: the full session in data/ (with the MATLAB result in
results/) if present, else the published excerpt in samples/ (scripts/make_samples.py)."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SAMPLE_OF = {
    "20261123_140216_ssbdet_v2": "ssbdet_pci165_150cuts",
    "20261123_112009_n1_2131": "n1_2131_20M_30bursts",
}


def find_session(name: str):
    """(session folder, expected-results CSV) or (None, None) if neither is there."""
    full, res = ROOT / "data" / name, ROOT / "results" / f"{name}_ssb.csv"
    if full.exists() and res.exists() and not os.environ.get("SNIFFER_SAMPLES_ONLY"):  # =1: as a fresh clone
        return full, res
    if name in SAMPLE_OF:
        s = ROOT / "samples" / SAMPLE_OF[name]
        if s.exists():
            return s, s / "expected_ssb.csv"
    return None, None
