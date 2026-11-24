"""Headless SSB monitor: decode every SSB from the board (or a replay), log it, raise alerts.

  .venv_local/bin/python -m capture.ssb_monitor --port /dev/ttyACM0 --freq 2131 --ssb 2125.45 \\
        --jsonl results/monitor.jsonl --csv results/monitor.csv --expect-pci 165 166
  .venv_local/bin/python -m capture.ssb_monitor --replay samples/ssbdet_pci165_150cuts --speed 0 --stdout
  .venv_local/bin/python -m capture.ssb_monitor ... --udp 127.0.0.1:5005      # JSON feed for other tools

Same board handling as the live view (capture/ssb_live.py): reconnects after board resets,
the chip's own cell search (firmware ssb3), decimated cuts, --record. Every decoded SSB becomes
one JSON object (type "ssb"); alerts and board events are objects of type "alert" / "event"
and are also printed to stderr:

  new_pci        a PCI decodes for the first time (alert if --expect-pci is given and it is
                 not in the list, else an event)
  mib_change     the MIB of a PCI changes (fields other than the SFN)
  sfn_jump       the SFN does not advance with the board's clock (more than 1 frame off)
  timing_jump    the SFN-anchored frame start leaves the running fit by more than --timing-alert-us
  cfo_outlier    the CFO of a decoded SSB is more than --cfo-alert-hz from the cell's median
  board          the board was lost and reconnected (with its reset reason, if the firmware tells)

Receive only. Alerts are hints, not proof: weak SSBs give outliers, and a board reset restarts
its clock (the timing state is reset then).
"""
from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import os
import queue
import socket
import sys
import threading
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .ssb_live import (MAX_BACKLOG, LiveSource, ReplaySource, _worker_init, add_source_args, check_source_args,
                       decode_job)

MIB_KEYS = ("scs_common_khz", "kssb", "dmrs_typeA_pos", "pdcch_config_sib1", "cell_barred", "intra_freq_reselection")
CSV_FIELDS = ("host_time", "t_pss_us", "pci", "ok", "sfn", "hrf", "issb", "snr_db", "rsrp_dbfs", "cfo_hz",
              "t_frame_us", "timing_resid_us", "pci_chip", "mode")


class CellState:
    """Per-PCI history for the alerts."""

    def __init__(self):
        self.mib = None
        self.cfo = deque(maxlen=50)
        self.tf = deque(maxlen=60)         # (t_frame_us, sfn)
        self.last = None                   # (t_pss_us, sfn)

    def check(self, r, a) -> list[tuple[str, str]]:
        out = []
        m = r["mib"]
        sig = {k: m[k] for k in MIB_KEYS}
        if self.mib is not None and sig != self.mib:
            changed = {k: (self.mib[k], sig[k]) for k in MIB_KEYS if self.mib[k] != sig[k]}
            out.append(("mib_change", f"PCI {r['pci']}: MIB changed {changed}"))
        self.mib = sig
        if len(self.cfo) >= 10 and abs(r["cfo"] - np.median(self.cfo)) > a.cfo_alert_hz:
            out.append(("cfo_outlier", f"PCI {r['pci']}: CFO {r['cfo']:+.0f} Hz, median {np.median(self.cfo):+.0f} Hz"))
        self.cfo.append(r["cfo"])
        tp, sfn = r.get("t_pss_us"), m["sfn"]
        if tp is not None and self.last is not None:
            frames = (tp - self.last[0]) / 1e4                       # 10 ms frames on the board's clock
            exp = (self.last[1] + round(frames)) % 1024
            if abs(frames - round(frames)) < 0.3 and min((sfn - exp) % 1024, (exp - sfn) % 1024) > 1:
                out.append(("sfn_jump", f"PCI {r['pci']}: SFN {sfn}, expected {exp} from the time"))
        if tp is not None:
            self.last = (tp, sfn)
        r["timing_resid_us"] = None
        if r.get("t_frame_us") is not None:
            if len(self.tf) >= 8:
                tf = np.array([x[0] for x in self.tf]); s = np.array([x[1] for x in self.tf])
                t0 = tf - s * 1e4
                t0 = t0 + np.round((t0[0] - t0) / 10.24e6) * 10.24e6
                p = np.polyfit(tf - tf[0], t0, 1)
                now0 = r["t_frame_us"] - sfn * 1e4
                now0 += np.round((t0[0] - now0) / 10.24e6) * 10.24e6
                resid = now0 - np.polyval(p, r["t_frame_us"] - tf[0])
                r["timing_resid_us"] = float(resid)
                if abs(resid) > a.timing_alert_us:
                    out.append(("timing_jump", f"PCI {r['pci']}: frame start {resid:+.2f} us off the fit"))
            self.tf.append((r["t_frame_us"], sfn))
        return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_source_args(p, label="monitor")
    p.add_argument("--jsonl", help="append every SSB, event and alert as JSON lines to this file")
    p.add_argument("--csv", help="append one CSV row per decoded SSB")
    p.add_argument("--udp", help="send every JSON object as a UDP datagram to host:port")
    p.add_argument("--stdout", action="store_true", help="print every JSON object to stdout")
    p.add_argument("--expect-pci", type=int, nargs="*", default=None, help="known PCIs; others raise an alert")
    p.add_argument("--cfo-alert-hz", type=float, default=2000.0)
    p.add_argument("--timing-alert-us", type=float, default=2.0)
    p.add_argument("--duration", type=float, default=None, help="seconds, then stop")
    a = p.parse_args(argv)
    check_source_args(p, a)

    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    pool = ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context("spawn"), initializer=_worker_init)
    events: queue.Queue = queue.Queue()
    results: queue.Queue = queue.Queue()
    lock, inflight = threading.Lock(), [0]

    def busy():
        return inflight[0] + events.qsize() >= MAX_BACKLOG

    src = ReplaySource(a, events, pool, busy) if a.replay else LiveSource(a, events, pool)
    jf = open(a.jsonl, "a") if a.jsonl else None
    cf = None
    if a.csv:
        new = not os.path.exists(a.csv) or os.path.getsize(a.csv) == 0
        cf = open(a.csv, "a", newline="")
        cw = csv.DictWriter(cf, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if new:
            cw.writeheader()
    udp = None
    if a.udp:
        host, port = a.udp.rsplit(":", 1)
        udp = (socket.socket(socket.AF_INET, socket.SOCK_DGRAM), (host, int(port)))

    def emit(obj):
        line = json.dumps(obj)
        if jf:
            jf.write(line + "\n"); jf.flush()
        if a.stdout:
            print(line, flush=True)
        if udp:
            try:
                udp[0].sendto(line.encode(), udp[1])
            except OSError:
                pass
        if obj["type"] in ("alert", "event"):
            print(f"[{time.strftime('%H:%M:%S')}] {obj['type'].upper()} {obj.get('kind', '')}: {obj['text']}",
                  file=sys.stderr, flush=True)

    def on_done(fut):
        try:
            results.put(fut.result())
        except Exception as e:
            results.put(e)
        with lock:
            inflight[0] -= 1

    cells: dict[int, CellState] = {}
    n = {"frames": 0, "decoded": 0, "mib": 0, "skipped": 0}
    last_chip = None
    src.start()
    t_end = time.monotonic() + a.duration if a.duration else None
    stopped = False
    try:
        while True:
            if t_end and time.monotonic() > t_end and not stopped:
                src.stop_flag.set(); stopped = True
            try:
                ev = events.get(timeout=0.05)
            except queue.Empty:
                ev = None
            if ev is not None:
                k = ev["kind"]
                if k == "frame":
                    n["frames"] += 1
                    tc = ev.get("t_pss_chip_us")
                    if tc is not None and last_chip is not None and tc < last_chip - 1e6:
                        for c in cells.values():           # the board's clock restarted (reset)
                            c.tf.clear(); c.last = None
                        emit({"type": "event", "kind": "board", "text": "board clock restarted: timing state reset"})
                    if tc is not None:
                        last_chip = tc
                    if ev.get("iq") is not None:
                        if not a.replay and inflight[0] >= MAX_BACKLOG:
                            n["skipped"] += 1
                        else:
                            job = {key: ev.get(key) for key in ("iq", "fs", "t_pss", "nid2", "cfo", "f_ssb", "f_lo",
                                                                "t_first_us", "seq", "mode", "stat", "R", "cfo_ddc",
                                                                "pci_chip", "scale")}
                            with lock:
                                inflight[0] += 1
                            pool.submit(decode_job, job).add_done_callback(on_done)
                elif k == "reconnect":
                    emit({"type": "alert", "kind": "board", "text": "board lost; reconnecting"})
                elif k == "log":
                    if "start reason" in ev["text"] or "board:" in ev["text"]:
                        emit({"type": "event", "kind": "board", "text": ev["text"]})
                elif k == "error":
                    emit({"type": "alert", "kind": "error", "text": ev["text"]})
                elif k == "stopped":
                    stopped = True
            try:
                while True:
                    item = results.get_nowait()
                    if isinstance(item, Exception):
                        emit({"type": "alert", "kind": "decode_error", "text": f"{type(item).__name__}: {item}"})
                        continue
                    for r in item:
                        if r.get("pci") is None:
                            continue
                        n["decoded"] += 1
                        rec = {"type": "ssb", "host_time": time.time(), "pci": r["pci"], "ok": bool(r["ok"]),
                               "issb": r.get("issb"), "snr_db": r.get("snr_db"), "rsrp_dbfs": r.get("rsrp_dbfs"),
                               "cfo_hz": r.get("cfo"), "t_pss_us": r.get("t_pss_us"), "t_frame_us": r.get("t_frame_us"),
                               "pci_chip": r.get("pci_chip"), "mode": r.get("mode"), "sss_stat": r.get("sss_stat")}
                        if r["ok"]:
                            n["mib"] += 1
                            m = r["mib"]
                            rec.update(sfn=m["sfn"], hrf=m["hrf"], **{k2: m[k2] for k2 in MIB_KEYS})
                            new_cell = r["pci"] not in cells
                            c = cells.setdefault(r["pci"], CellState())
                            if new_cell:
                                unexpected = a.expect_pci is not None and r["pci"] not in a.expect_pci
                                emit({"type": "alert" if unexpected else "event", "kind": "new_pci",
                                      "text": f"PCI {r['pci']} decoded" + (" (not in --expect-pci)" if unexpected else ""),
                                      "pci": r["pci"]})
                            for kind, text in c.check(r, a):
                                emit({"type": "alert", "kind": kind, "text": text, "pci": r["pci"]})
                            rec["timing_resid_us"] = r.get("timing_resid_us")
                        emit(rec)
                        if cf:
                            cw.writerow(rec); cf.flush()
            except queue.Empty:
                pass
            if stopped and not src.is_alive() and inflight[0] == 0 and events.empty() and results.empty():
                break
    except KeyboardInterrupt:
        src.stop_flag.set()
        src.join(timeout=8)
    finally:
        pool.shutdown(wait=False, cancel_futures=True)
        for fh in (jf, cf):
            if fh:
                fh.close()
    print(f"ssb_monitor: {n['frames']} frames, {n['decoded']} decoded, {n['mib']} MIBs, {n['skipped']} skipped; "
          f"cells {sorted(cells)}", file=sys.stderr)


if __name__ == "__main__":
    main()
