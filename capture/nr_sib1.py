"""SIB1 (TS 38.331 BCCH-DL-SCH-Message, unaligned PER) decoded with pycrate, and what it says in words.

  .venv_local/bin/python -m capture.nr_sib1 results/sib1_<session>_<PCI>.json [--asn1]

The input is the JSON that matlab/sib1_offline.m writes (decoded SIB1 transport blocks as hex).
pycrate (LGPL, TS 38.331 ASN.1 included) is an optional dependency: pip install pycrate.

Frequencies (TS 38.211 4.4.4.2, 7.4.3.1; TS 38.104 5.4.2.1):
- the SSB's subcarrier 0 is SS_REF - 120 x 15 kHz; the common resource block holding it starts
  k_SSB x 15 kHz lower (subCarrierSpacingCommon 15 kHz);
- downlink point A lies offsetToPointA x 180 kHz below that (FrequencyInfoDL-SIB);
- uplink point A is given as an NR-ARFCN (absoluteFrequencyPointA);
- each carrier's resource blocks start offsetToCarrier RBs above point A, carrierBandwidth RBs wide.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SCS_HZ = {"kHz15": 15e3, "kHz30": 30e3, "kHz60": 60e3, "kHz120": 120e3}


def nr_arfcn_to_hz(n: int) -> float:
    """NR-ARFCN -> frequency (TS 38.104 Table 5.4.2.1-1)."""
    if n < 600000:
        return n * 5e3
    if n < 2016667:
        return 3000e6 + (n - 600000) * 15e3
    return 24250.08e6 + (n - 2016667) * 60e3


def sliv(v: int) -> tuple[int, int]:
    """startSymbolAndLength -> (S, L), TS 38.214 5.1.2.1: SLIV = 14 (L-1) + S if L-1 <= 7, else
    14 (14-L+1) + (14-1-S), with S + L <= 14."""
    a, b = divmod(v, 14)
    if a + 1 + b <= 14 and a <= 7:
        return b, a + 1
    return 13 - b, 15 - a


def matlab_options(sib1: dict, coreset0_offset_rb: int | None = None) -> str:
    """The options matlab/sib1_offline.m needs for SI messages other than SIB1."""
    dl = sib1["servingCellConfigCommon"]["downlinkConfigCommon"]
    opts = [f"'OffsetToPointA', {dl['frequencyInfoDL']['offsetToPointA']}"]
    pc = dl["initialDownlinkBWP"].get("pdsch-ConfigCommon")
    if pc and pc[0] == "setup" and "pdsch-TimeDomainAllocationList" in pc[1]:
        rows = []
        for e in pc[1]["pdsch-TimeDomainAllocationList"]:
            S, L = sliv(e["startSymbolAndLength"])
            rows.append(f"{1 if e['mappingType'] == 'typeA' else 0} {S} {L}")
        opts.append(f"'TimeAllocationList', [{'; '.join(rows)}]")
    return ", ".join(opts)


WINDOW_SLOTS = {f"s{n}": n for n in (5, 10, 20, 40, 80, 160, 320, 640, 1280)}


def si_plan(sib1: dict, slots_per_frame: int = 10) -> dict | None:
    """Where the SI messages other than SIB1 are sent, for a capture (capture/sib1_capture.py) and for
    matlab/sib1_offline.m; None if SIB1 schedules none.

    TS 38.331 5.2.2.3.2: SI message n (from 0, in schedulingInfoList order) has a window of w slots
    (si-WindowLength) that starts x = n w slots after the start of every frame with SFN mod T = 0
    (T = si-Periodicity in frames). Its PDCCH is in searchSpaceOtherSystemInformation: 0 (or absent)
    means SIB1's Type0 occasions, another id one of the commonSearchSpaceList (slots from its
    monitoringSlotPeriodicityAndOffset; all slots if the period does not divide a frame).
    Returns {"messages": [{"sibs", "T", "frames"}], "sfn_mod", "sfn_set" (frames to capture, SFN mod
    sfn_mod), "slots" (None: SIB1's occasions), "window_slots", "search_space", "matlab" (options)}."""
    si = sib1.get("si-SchedulingInfo")
    if not si or not si.get("schedulingInfoList"):
        return None
    w = WINDOW_SLOTS[si["si-WindowLength"]]
    msgs = []
    for n, s in enumerate(si["schedulingInfoList"]):
        x = n * w
        msgs.append({"sibs": [m["type"].replace("sibType", "SIB") for m in s["sib-MappingInfo"]],
                     "T": int(s["si-Periodicity"].replace("rf", "")), "status": s.get("si-BroadcastStatus"),
                     "frames": sorted({(x + k) // slots_per_frame for k in range(w)})})
    sfn_mod = math.lcm(*(m["T"] for m in msgs))
    sfn_set = sorted({(f + j * m["T"]) % sfn_mod for m in msgs for f in m["frames"] for j in range(sfn_mod // m["T"])})
    bwp = sib1["servingCellConfigCommon"]["downlinkConfigCommon"]["initialDownlinkBWP"]
    pc = bwp.get("pdcch-ConfigCommon")
    pc = pc[1] if pc and pc[0] == "setup" else {}
    ss_id = pc.get("searchSpaceOtherSystemInformation", 0)
    slots = None
    if ss_id:
        ss = next((s for s in pc.get("commonSearchSpaceList", []) if s["searchSpaceId"] == ss_id), None)
        slots = list(range(slots_per_frame))
        if ss and "monitoringSlotPeriodicityAndOffset" in ss:
            per, off = ss["monitoringSlotPeriodicityAndOffset"]
            p = int(per.replace("sl", ""))
            if slots_per_frame % p == 0:
                slots = [k for k in range(slots_per_frame) if (k - off) % p == 0]
    opts = matlab_options(sib1)
    if slots is not None:
        opts += (f", 'Slots', [{' '.join(map(str, slots))}], 'FrameMod', {sfn_mod}, "
                 f"'FrameSet', [{' '.join(map(str, sfn_set))}]")
    return {"messages": msgs, "sfn_mod": sfn_mod, "sfn_set": sfn_set, "slots": slots, "window_slots": w,
            "search_space": ss_id, "matlab": opts}


def decode(payload: bytes) -> dict:
    """BCCH-DL-SCH-Message -> {"text": ASN.1 value notation, "sib1": the SIB1 value (dict) or None,
    "sibs": [(type, value)] of a SystemInformation message, or []}."""
    try:
        from pycrate_asn1dir import RRCNR
    except ImportError as e:
        raise ImportError("SIB1 decoding needs pycrate: .venv_local/bin/pip install pycrate") from e
    m = RRCNR.NR_RRC_Definitions.BCCH_DL_SCH_Message
    m.from_uper(payload)
    v = m.get_val()["message"]
    sib1 = v[1][1] if v[0] == "c1" and v[1][0] == "systemInformationBlockType1" else None
    sibs = []
    if v[0] == "c1" and v[1][0] == "systemInformation":
        ce = v[1][1]["criticalExtensions"]
        if ce[0] == "systemInformation":
            sibs = list(ce[1]["sib-TypeAndInfo"])
    return {"text": m.to_asn1(), "sib1": sib1, "sibs": sibs}


def _carrier(point_a_hz: float, scs_list: list) -> list[tuple]:
    out = []
    for c in scs_list:
        rb = 12 * SCS_HZ[c["subcarrierSpacing"]]
        lo = point_a_hz + c["offsetToCarrier"] * rb
        out.append((c["subcarrierSpacing"], c["carrierBandwidth"], lo, lo + c["carrierBandwidth"] * rb))
    return out


def summary(sib1: dict, f_ssb_hz: float | None = None, kssb: int | None = None) -> list[tuple[str, str]]:
    """(label, text) rows: who runs the cell, where its carriers are, what else it broadcasts."""
    rows = []
    for k, p in enumerate(sib1["cellAccessRelatedInfo"]["plmn-IdentityInfoList"]):
        plmns = ", ".join("".join(map(str, x.get("mcc", []))) + "-" + "".join(map(str, x["mnc"]))
                          for x in p["plmn-IdentityList"])
        rows.append((f"network (PLMN) {k + 1}", plmns))
        if "trackingAreaCode" in p:
            tac = p["trackingAreaCode"][0]
            rows.append(("tracking area code", f"{tac} (0x{tac:06X})"))
        nci = p["cellIdentity"][0]
        rows.append(("cell identity (NCI, 36 bit)", f"0x{nci:09X}"))
        rows.append(("reserved for operator use", p["cellReservedForOperatorUse"]))
    cs = sib1.get("cellSelectionInfo", {})
    if "q-RxLevMin" in cs:
        rows.append(("min. RX level (q-RxLevMin)", f"{2 * cs['q-RxLevMin']} dBm"))
    scc = sib1.get("servingCellConfigCommon")
    if scc:
        dl = scc["downlinkConfigCommon"]["frequencyInfoDL"]
        rows.append(("bands", ", ".join(f"n{b['freqBandIndicatorNR']}" for b in dl["frequencyBandList"])))
        if f_ssb_hz is not None and kssb is not None:
            crb = f_ssb_hz - 120 * 15e3 - kssb * 15e3
            pa = crb - dl["offsetToPointA"] * 180e3
            rows.append(("downlink point A", f"{pa / 1e6:.3f} MHz (offsetToPointA {dl['offsetToPointA']} RB)"))
            for scs, n, lo, hi in _carrier(pa, dl["scs-SpecificCarrierList"]):
                rows.append(("downlink carrier", f"{n} RB at {scs}: {lo / 1e6:.3f}-{hi / 1e6:.3f} MHz "
                                                 f"(centre {(lo + hi) / 2e6:.3f} MHz, {(hi - lo) / 1e6:.2f} MHz of RBs)"))
        ul = scc.get("uplinkConfigCommon", {}).get("frequencyInfoUL")
        if ul and "absoluteFrequencyPointA" in ul:
            pa = nr_arfcn_to_hz(ul["absoluteFrequencyPointA"])
            for scs, n, lo, hi in _carrier(pa, ul["scs-SpecificCarrierList"]):
                rows.append(("uplink carrier", f"{n} RB at {scs}: {lo / 1e6:.3f}-{hi / 1e6:.3f} MHz "
                                               f"(centre {(lo + hi) / 2e6:.3f} MHz)"))
            if "p-Max" in ul:
                rows.append(("max. UE power (p-Max)", f"{ul['p-Max']} dBm"))
        pos = scc.get("ssb-PositionsInBurst", {})
        if "inOneGroup" in pos:
            v, nb = pos["inOneGroup"]
            idx = [i for i in range(nb) if (v >> (nb - 1 - i)) & 1]
            rows.append(("SSB indices sent", ", ".join(map(str, idx))))
        if "ssb-PeriodicityServingCell" in scc:
            rows.append(("SSB period", scc["ssb-PeriodicityServingCell"].replace("ms", "") + " ms"))
        if "ss-PBCH-BlockPower" in scc:
            rows.append(("SSB power (EPRE)", f"{scc['ss-PBCH-BlockPower']} dBm per RE"))
        tdd = scc.get("tdd-UL-DL-ConfigurationCommon")
        rows.append(("duplex", "TDD" if tdd else "FDD (no TDD pattern)"))
    si = sib1.get("si-SchedulingInfo")
    if si:
        for s in si["schedulingInfoList"]:
            sibs = ", ".join(x["type"].replace("sibType", "SIB") for x in s["sib-MappingInfo"])
            rows.append(("other SIBs", f"{sibs} every {s['si-Periodicity'].replace('rf', '')} frames "
                                       f"({s['si-BroadcastStatus']})"))
    for k in ("ims-EmergencySupport", "eCallOverIMS-Support"):
        if k in sib1:
            rows.append((k, str(sib1[k])))
    return rows


# LTE downlink EARFCN -> MHz (TS 36.101 Table 5.7.3-1): band: (F_DL_low MHz, N_Offs-DL, last EARFCN), for the
# FDD/TDD bands used in Europe; other EARFCNs are shown as numbers only
EUTRA_DL = {1: (2110.0, 0, 599), 3: (1805.0, 1200, 1949), 7: (2620.0, 2750, 3449), 8: (925.0, 3450, 3799),
            20: (791.0, 6150, 6449), 28: (758.0, 9210, 9659), 32: (1452.0, 9920, 10359), 38: (2570.0, 37750, 38249),
            40: (2300.0, 38650, 39649)}
MBW_MHZ = {"mbw6": 1.4, "mbw15": 3, "mbw25": 5, "mbw50": 10, "mbw75": 15, "mbw100": 20}


def earfcn_to_mhz(n: int) -> tuple[int, float] | None:
    for band, (f0, off, last) in EUTRA_DL.items():
        if off <= n <= last:
            return band, f0 + 0.1 * (n - off)
    return None


def _neigh(lst) -> str:
    return ", ".join(str(c["physCellId"]) + (f" ({c['q-OffsetCell']})" if c.get("q-OffsetCell", "dB0") != "dB0" else "")
                     for c in lst)


def summary_sib(typ: str, val) -> list[tuple[str, str]]:
    """SIB2-SIB5 in words: when a phone reselects, and which other cells and carriers exist."""
    rows = []
    if typ == "sib2":
        c = val.get("cellReselectionInfoCommon", {})
        if "q-Hyst" in c:
            rows.append(("reselection hysteresis", c["q-Hyst"].replace("dB", "") + " dB"))
        f = val.get("cellReselectionServingFreqInfo", {})
        if "cellReselectionPriority" in f:
            rows.append(("priority of this frequency", str(f["cellReselectionPriority"])))
        i = val.get("intraFreqCellReselectionInfo", {})
        if "s-IntraSearchP" in i:
            rows.append(("measure neighbours below", f"{2 * i['s-IntraSearchP']} dB over q-RxLevMin (s-IntraSearchP)"))
        sm = i.get("smtc", {})
        if sm:
            per, off = sm["periodicityAndOffset"]
            rows.append(("SSB measurement window", f"{sm['duration'].replace('sf', '')} ms every {per.replace('sf', '')} ms, "
                                                   f"offset {off} ms (SMTC)"))
    elif typ == "sib3":
        lst = val.get("intraFreqNeighCellList", [])
        rows.append((f"neighbour cells, same freq. ({len(lst)})", _neigh(lst) or "-"))
        if val.get("intraFreqBlackCellList"):
            rows.append(("excluded cells", str(val["intraFreqBlackCellList"])))
    elif typ == "sib4":
        for c in val.get("interFreqCarrierFreqList", []):
            f = nr_arfcn_to_hz(c["dl-CarrierFreq"]) / 1e6
            bands = "/".join(f"n{b['freqBandIndicatorNR']}" for b in c.get("frequencyBandList", []))
            rows.append((f"NR carrier {bands}", f"{f:.2f} MHz (NR-ARFCN {c['dl-CarrierFreq']}), SSB "
                                                f"{c.get('ssbSubcarrierSpacing', '?').replace('kHz', '')} kHz, priority "
                                                f"{c.get('cellReselectionPriority', '-')}; cells "
                                                f"{_neigh(c.get('interFreqNeighCellList', [])) or '-'}"))
    elif typ == "sib5":
        for c in val.get("carrierFreqListEUTRA", []):
            bf = earfcn_to_mhz(c["carrierFreq"])
            where = f"band {bf[0]}, {bf[1]:.1f} MHz" if bf else "band not in the table"
            bw = MBW_MHZ.get(c.get("allowedMeasBandwidth"), "?")
            rows.append(("LTE carrier", f"EARFCN {c['carrierFreq']} ({where}), {bw} MHz, priority "
                                        f"{c.get('cellReselectionPriority', '-')}"))
    else:
        rows.append((typ, "decoded; no summary (use --asn1)"))
    return rows


_DECODED: dict[str, dict] = {}                     # hex -> decode() (pycrate is slow to load and to run)


def load_decoded(pci: int, ssb_mhz: float, folder: Path | None = None) -> dict:
    """Everything decoded so far for one cell: the results/sib1_*_<pci>.json of matlab/sib1_offline.m with
    the same SSB frequency, newest first. Returns {"sib1": value or None, "kssb", "sib1_rows": summary,
    "sib1_from": (session, copies), "sibs": {type: (value, rows, session)}, "plan": si_plan or None,
    "files": [...], "error": text or None}. The newest decode of each message wins."""
    folder = folder or ROOT / "results"
    out = {"sib1": None, "kssb": None, "sib1_rows": [], "sib1_from": None, "sibs": {}, "plan": None, "files": [],
           "error": None}
    files = []
    for f in folder.glob(f"sib1_*_{pci}.json"):
        try:
            j = json.loads(f.read_text())
        except (OSError, ValueError):
            continue
        if j.get("pci") == pci and abs(float(j.get("ssb_mhz", 0)) - ssb_mhz) < 0.01:
            files.append((j.get("session", f.stem), f, j))
    files.sort(key=lambda t: t[0], reverse=True)
    out["files"] = [str(f) for _, f, _ in files]
    for session, f, j in files:
        items = j["sib1"] if isinstance(j["sib1"], list) else ([j["sib1"]] if j.get("sib1") else [])
        copies: dict[str, int] = {}
        for it in items:
            copies[it["hex"]] = copies.get(it["hex"], 0) + 1
        for hx, n in copies.items():
            try:
                d = _DECODED.get(hx) or _DECODED.setdefault(hx, decode(bytes.fromhex(hx)))
            except ImportError as e:
                out["error"] = str(e)
                return out
            except Exception as e:  # a transport block that passed the CRC but is no valid message
                out["error"] = f"{session}: {type(e).__name__}: {e}"
                continue
            if d["sib1"] is not None and out["sib1"] is None:
                out.update(sib1=d["sib1"], kssb=j.get("kssb"), sib1_from=(session, n),
                           sib1_rows=summary(d["sib1"], ssb_mhz * 1e6, j.get("kssb")))
            for typ, val in d["sibs"]:
                if typ not in out["sibs"]:
                    out["sibs"][typ] = (val, summary_sib(typ, val), session)
    if out["sib1"] is not None:
        out["plan"] = si_plan(out["sib1"])
    return out


def neighbour_pcis(info: dict) -> list[int]:
    """Same-frequency neighbour PCIs from a decoded SIB3 (load_decoded), in list order."""
    sib3 = info.get("sibs", {}).get("sib3")
    return [c["physCellId"] for c in sib3[0].get("intraFreqNeighCellList", [])] if sib3 else []


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("json", help="results/sib1_<session>_<PCI>.json from matlab/sib1_offline.m")
    p.add_argument("--kssb", type=int, default=None, help="k_SSB from the MIB (default: from the JSON)")
    p.add_argument("--asn1", action="store_true", help="also print the whole message in ASN.1 value notation")
    p.add_argument("--matlab", action="store_true", help="print the sib1_offline options for the other SI messages")
    a = p.parse_args(argv)
    j = json.load(open(a.json))
    items = j["sib1"] if isinstance(j["sib1"], list) else [j["sib1"]]
    kssb = a.kssb if a.kssb is not None else j.get("kssb")
    seen = {}
    for it in items:
        seen.setdefault(it["hex"], []).append(f"SFN {it['sfn']} slot {it['slot']}")
    if not seen:
        print("no decoded SIB1 in", a.json)
        return 1
    for hx, where in seen.items():
        d = decode(bytes.fromhex(hx))
        what = "SIB1" if d["sib1"] is not None else "SI message " + "+".join(t.upper() for t, _ in d["sibs"])
        print(f"PCI {j['pci']}, SSB {j['ssb_mhz']} MHz: {what}, {len(hx) // 2} bytes, decoded {len(where)}x "
              f"({', '.join(sorted(set(where)))})")
        if d["sib1"] is not None:
            for label, text in summary(d["sib1"], j["ssb_mhz"] * 1e6, kssb):
                print(f"  {label:28s} {text}")
            if a.matlab:
                print("  sib1_offline options for SI messages:", matlab_options(d["sib1"]))
        for typ, val in d["sibs"]:
            print(f"  -- {typ} --")
            for label, text in summary_sib(typ, val):
                print(f"  {label:28s} {text}")
        if a.asn1:
            print(d["text"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
