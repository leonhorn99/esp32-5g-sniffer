"""What a decoded MIB says, in words: the fields of nr_ssb.parse_mib, where the cell's SIB1 control
channel is (CORESET#0 and its search space 0), and short explanations for the live view.

Tables from 3GPP TS 38.213 V17.14.0 (ETSI TS 138 213, 2026-09), clause 13, copied from the ETSI PDF
and checked row by row on 24 Nov 2026:
  Table 13-1   CORESET for Type0-PDCCH, {SS/PBCH block, PDCCH} SCS {15, 15} kHz, bands with minimum
               channel bandwidth 5 or 10 MHz (n1 and n7 are such bands)
  Table 13-11  Type0-PDCCH monitoring occasions, multiplexing pattern 1, FR1
Other SCS pairs (Table 13-2 ff.) are not included: the live decoder handles 15 kHz SSBs only, and a
cell with subCarrierSpacingCommon 30 kHz is reported as such.
"""
from __future__ import annotations

from fractions import Fraction

SCS_KHZ = 15
# Table 13-1: index -> (multiplexing pattern, N_RB, N_symb, offset in RBs); 15 is reserved
CORESET0_15_15 = {
    0: (1, 24, 2, 0), 1: (1, 24, 2, 2), 2: (1, 24, 2, 4), 3: (1, 24, 3, 0), 4: (1, 24, 3, 2), 5: (1, 24, 3, 4),
    6: (1, 48, 1, 12), 7: (1, 48, 1, 16), 8: (1, 48, 2, 12), 9: (1, 48, 2, 16), 10: (1, 48, 3, 12),
    11: (1, 48, 3, 16), 12: (1, 96, 1, 38), 13: (1, 96, 2, 38), 14: (1, 96, 3, 38),
}
# Table 13-11: index -> (O, search space sets per slot, M, first symbol: 0, 1, 2 or "alt" =
# {0 if i even, N_symb^CORESET if i odd})
SS0_PATTERN1_FR1 = {
    0: (0, 1, Fraction(1), 0), 1: (0, 2, Fraction(1, 2), "alt"), 2: (2, 1, Fraction(1), 0),
    3: (2, 2, Fraction(1, 2), "alt"), 4: (5, 1, Fraction(1), 0), 5: (5, 2, Fraction(1, 2), "alt"),
    6: (7, 1, Fraction(1), 0), 7: (7, 2, Fraction(1, 2), "alt"), 8: (0, 1, Fraction(2), 0),
    9: (5, 1, Fraction(2), 0), 10: (0, 1, Fraction(1), 1), 11: (0, 1, Fraction(1), 2),
    12: (2, 1, Fraction(1), 1), 13: (2, 1, Fraction(1), 2), 14: (5, 1, Fraction(1), 1),
    15: (5, 1, Fraction(1), 2),
}

# (key, label, explanation) of the rows of the live view's MIB tab, in order
FIELDS = (
    ("n", "MIBs decoded",
     "How many SSBs of this cell passed the PBCH CRC since the stream (or channel) started. Only a passed CRC "
     "proves the cell is real, not a noise peak with a random PCI."),
    ("content", "MIB content",
     "Whether any field other than SFN and half frame changed between this cell's MIBs. A change means the "
     "cell was reconfigured."),
    ("sfn", "SFN",
     "System frame number, 0-1023, +1 every 10 ms (6 bits in the MIB, 4 more in the PBCH payload). It places "
     "each SSB on the cell's frame grid; the frame-timing plot is built on it."),
    ("hrf", "half frame",
     "0: the SSB was in the first 5 ms of the 10 ms frame, 1: in the second."),
    ("issb", "SSB index (beam)",
     "Which SSB of the cell's burst set, i.e. which beam. It comes from the PBCH DM-RS sequence (with L_max 4: "
     "2 bits, plus the half frame). One index means we see one beam, or the cell sends only one."),
    ("scs", "subCarrierSpacingCommon",
     "Subcarrier spacing of SIB1, paging and the random-access messages, i.e. of initial access."),
    ("kssb", "k_SSB (subcarrier offset)",
     "Offset in 15 kHz subcarriers between the SSB and the carrier's common resource-block grid "
     "(ssb-SubcarrierOffset plus one bit from the PBCH payload). Below 24 a CORESET for SIB1 exists "
     "(TS 38.213 4.1)."),
    ("sib1", "SIB1 announced",
     "k_SSB < 24: the MIB points to the control channel of SIB1, so the cell broadcasts SIB1, as a cell that "
     "phones can attach to directly (standalone) does. Cells that only add capacity to LTE (non-standalone) "
     "often send none. The SA/NSA reading is an inference."),
    ("dmrs", "dmrs-TypeA-Position",
     "OFDM symbol (2 or 3) of the first DM-RS of the data channel (PDSCH) in each slot."),
    ("pdcch", "pdcch-ConfigSIB1",
     "8 bits: upper 4 = controlResourceSetZero (the CORESET#0 row of TS 38.213 Table 13-1), lower 4 = "
     "searchSpaceZero (when to look, Table 13-11)."),
    ("coreset0", "CORESET#0 size",
     "The control region that carries the PDCCH scheduling SIB1: resource blocks x OFDM symbols, from TS 38.213 "
     "Table 13-1 ({SSB, PDCCH} SCS {15, 15} kHz)."),
    ("coreset0_mhz", "CORESET#0 band (MHz)",
     "Computed from the MIB: its lowest RB lies the Table 13-1 offset (in RBs) below the common RB that holds "
     "the SSB's first subcarrier (SS_REF - 120 x 15 kHz - k_SSB x 15 kHz); width N_RB x 180 kHz."),
    ("ss0_slots", "SIB1 PDCCH slots",
     "TS 38.213 Table 13-11 (pattern 1, FR1) gives O and M. For SSB index i the PDCCH for SIB1 is in slots n0 "
     "and n0+1 (of 10 per frame at 15 kHz), n0 = (O + floor(i x M)) mod 10, of frames with SFN mod 2 = "
     "floor((O + floor(i x M)) / 10) mod 2. Shown for the most frequent SSB index."),
    ("ss0_sym", "SIB1 PDCCH symbol",
     "First OFDM symbol of CORESET#0 in those slots (Table 13-11)."),
    ("barred", "cellBarred",
     "barred: phones must not camp on this cell."),
    ("intra", "intraFreqReselection",
     "Whether phones, when this cell is barred, may reselect another cell on the same frequency."),
)
HELP = {k: h for k, _, h in FIELDS}


def content_key(mib: dict) -> tuple:
    """The MIB fields that stay fixed while the cell is configured the same (all but SFN, half frame)."""
    return tuple(sorted((k, v) for k, v in mib.items() if k not in ("sfn", "hrf", "path_metric", "bits")))


def coreset0(mib: dict) -> dict | None:
    """Table 13-1 row of the MIB, or None (SCS pair not included, reserved index, or no CORESET)."""
    if mib["kssb"] >= 24 or mib["scs_common_khz"] != SCS_KHZ:
        return None
    row = CORESET0_15_15.get(mib["pdcch_config_sib1"] >> 4)
    if row is None:
        return None
    pattern, n_rb, n_sym, offset = row
    return {"pattern": pattern, "n_rb": n_rb, "n_sym": n_sym, "offset_rb": offset}


def coreset0_band_hz(mib: dict, f_ssb_hz: float) -> tuple[float, float] | None:
    """Lowest and highest frequency (Hz) of CORESET#0 for an SSB at SS_REF f_ssb_hz (15 kHz SSB, 15 kHz CORESET)."""
    c = coreset0(mib)
    if c is None or f_ssb_hz is None:
        return None
    scs = SCS_KHZ * 1e3
    crb0 = f_ssb_hz - 120 * scs - mib["kssb"] * scs       # subcarrier 0 of the common RB holding the SSB's first subcarrier
    lo = crb0 - c["offset_rb"] * 12 * scs
    return lo, lo + c["n_rb"] * 12 * scs


def sib1_occasions(mib: dict, issb: int) -> dict | None:
    """Type0-PDCCH monitoring for SSB index issb (Table 13-11, mu = 0): slots n0 and n0 + 1, frame parity,
    first symbol."""
    c = coreset0(mib)
    if c is None or c["pattern"] != 1:
        return None
    o, per_slot, m, first = SS0_PATTERN1_FR1[mib["pdcch_config_sib1"] & 15]
    n_slot = 10                                              # slots per frame at 15 kHz
    v = o * 1 + int(issb * m)                                # O x 2^mu + floor(i x M), mu = 0
    sym = (0 if issb % 2 == 0 else c["n_sym"]) if first == "alt" else first
    return {"O": o, "M": m, "per_slot": per_slot, "slots": (v % n_slot, v % n_slot + 1),
            "frame_parity": (v // n_slot) % 2, "first_symbol": sym}


def describe(mib: dict, f_ssb_hz: float | None, issb: int) -> dict:
    """Display strings for the rows of FIELDS that depend on one MIB (not n / content)."""
    p = mib["pdcch_config_sib1"]
    c = coreset0(mib)
    band = coreset0_band_hz(mib, f_ssb_hz)
    occ = sib1_occasions(mib, issb)
    if mib["kssb"] >= 24:
        cs = "none"
    elif c is None:
        cs = "not in table" if mib["scs_common_khz"] != SCS_KHZ else "reserved"
    else:
        cs = f"{c['n_rb']} RB × {c['n_sym']} sym"
    return {
        "sfn": str(mib["sfn"]), "hrf": str(mib["hrf"]),
        "scs": f"{mib['scs_common_khz']} kHz", "kssb": str(mib["kssb"]),
        "sib1": "yes" if mib["kssb"] < 24 else "no",
        "dmrs": f"pos{mib['dmrs_typeA_pos']}",
        "pdcch": f"{p} = {p >> 4} / {p & 15}",
        "coreset0": cs,
        "coreset0_mhz": f"{band[0] / 1e6:.2f}-{band[1] / 1e6:.2f}" if band else "-",
        "ss0_slots": f"{occ['slots'][0]}, {occ['slots'][1]} ({('even', 'odd')[occ['frame_parity']]} SFN)" if occ else "-",
        "ss0_sym": str(occ["first_symbol"]) if occ else "-",
        "barred": "barred" if mib["cell_barred"] else "notBarred",
        "intra": str(mib["intra_freq_reselection"]),
    }
