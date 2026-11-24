#!/usr/bin/env bash
# First light of the ESP32-C5 SSB sniffer, receive only, no input needed.
# Run after the board is flashed (firmware/fw.sh flash). From anywhere:
#
#   scripts/first_light.sh [PORT]          # default /dev/ttyACM0; COUNT=300 bursts per LO
#
# 1. checks that the board answers as ESP-SDR/ESPsoup (INFO, CAPS with CAPT)
# 2. Wi-Fi sanity capture at 2442 MHz (bursts must show packets: proves the receive path)
# 3. n1 over the air, 20 MS/s, 10 bit, automatic gain, at three LO positions:
#      2131 MHz  window 2123-2139: top of A1's block (2110-2135) + start of Magenta's
#      2142 MHz  window 2134-2150: Magenta's block (2135-2150), NR via DSS
#      2160 MHz  window 2152-2168: Drei's block (2150-2170)
# 4. runs matlab/detect_session.m (full GSCN search) on each n1 session
# Captures: data/<time>_<label>/, results: results/<session>_ssb.md, log: logs/first_light_<time>.log
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PORT=${1:-/dev/ttyACM0}
COUNT=${COUNT:-300}
PY=$ROOT/.venv_local/bin/python
STAMP=$(date +%Y%m%d_%H%M%S)
LOG=$ROOT/logs/first_light_$STAMP.log
mkdir -p "$ROOT/logs"
cd "$ROOT"
exec > >(tee -a "$LOG") 2>&1

echo "== 1. firmware check on $PORT"
if ! INFO=$(timeout 60 "$PY" -m capture.capture_bursts info --port "$PORT" 2>&1); then
    echo "$INFO" | tail -5
    echo "The board does not answer as ESP-SDR. Flash it first: firmware/fw.sh flash $PORT"
    exit 1
fi
echo "$INFO" | "$PY" -c 'import json,sys; d=json.load(sys.stdin); print(" ", d["info"], "|", d.get("fwopt"), "| CAPT" if "CAPT" in d["caps"] else "| no CAPT")'

capture() {   # capture <label> <count> <freq...>
    local label=$1 count=$2; shift 2
    "$PY" -m capture.capture_bursts capture --port "$PORT" --freq "$@" --rate 20 --bits 10 \
        --gain auto --count "$count" --label "$label" --note "first light $STAMP, board antenna, indoors"
    ls -td "$ROOT"/data/*_"$label" | head -1
}

echo "== 2. Wi-Fi sanity capture, 2442 MHz"
WIFI=$(capture wifi_2442 20 2442 | tail -1)
"$PY" - "$WIFI" <<'EOF'
import sys, numpy as np
from capture.session import load_session
meta, recs, iq = load_session(sys.argv[1])
p = [np.mean(np.abs(x)**2) for x in iq]
peak = [np.max(np.abs(x)**2) for x in iq]
ratio = np.median(np.array(peak)/np.array(p))
print(f"  {len(iq)} bursts, gain {meta['settings']['gain']}, median peak/mean power {10*np.log10(ratio):.1f} dB "
      f"({'packet bursts visible' if ratio > 20 else 'no clear packets: check antenna and gain'})")
EOF

SESSIONS=()
for LO in 2131 2142 2160; do
    echo "== 3. n1 capture, LO $LO MHz, $COUNT bursts"
    SESSIONS+=("$(capture "n1_$LO" "$COUNT" "$LO" | tail -1)")
done

echo "== 4. SSB detection (MATLAB, nice 19, 2 threads; about 1 s per burst)"
for S in "${SESSIONS[@]}"; do
    nice -n 19 matlab -nodisplay -nosplash -batch \
        "maxNumCompThreads(2); addpath('$ROOT/matlab'); detect_session('$S');" \
        > "$ROOT/logs/detect_$(basename "$S").log" 2>&1 < /dev/null || echo "  MATLAB failed, see logs/detect_$(basename "$S").log"
    tail -1 "$ROOT/logs/detect_$(basename "$S").log"
done

echo "== summaries"
for S in "${SESSIONS[@]}"; do
    R=$ROOT/results/$(basename "$S")_ssb.md
    [[ -f $R ]] && sed -n '1,30p' "$R" && echo
done
echo "log: $LOG"
