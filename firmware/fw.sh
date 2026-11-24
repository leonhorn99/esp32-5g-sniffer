#!/usr/bin/env bash
# Build and flash the ESPsoup / ESP-SDR firmware for the ESP32-C5 with the pinned ESP-IDF.
#
#   firmware/fw.sh build                    # build firmware/espsoup-firmware (branch ssb-sniffer)
#   firmware/fw.sh flash [PORT]             # flash the build (default PORT /dev/ttyACM0)
#   firmware/fw.sh erase [PORT]             # erase the whole flash (back it up first, see README)
#   firmware/fw.sh flash-prebuilt [PORT]    # flash an official ESPsoup 0.7.1 image put in firmware/prebuilt/
#   firmware/fw.sh info [PORT]              # ask the running firmware (INFO, CAPS, LIMITS, VERSION)
#   firmware/fw.sh env                      # print the ESP-IDF and build paths in use
#
# The firmware calls undocumented PHY functions from Espressif's binary libraries, so it
# is built with exactly the ESP-IDF commit its firmware-targets.json pins (as its CI does).
# That IDF lives in ~/.espressif/idf-<commit8>/esp-idf, next to any other ESP-IDF install.
# The build directory is outside the source tree (thousands of files): ~/build/esp32_5g_sniffer/.
# Overrides: ESPSDR_IDF_PATH, ESPSDR_BUILD_DIR.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(dirname "$HERE")
SRC=$HERE/espsoup-firmware
PREBUILT=$HERE/prebuilt/espsoup-0.7.1
PROFILE=esp32c5
PIN=$(python3 -c "import json,sys; p={x['id']:x for x in json.load(open(sys.argv[1]))['profiles']}; print(p[sys.argv[2]]['idf_ref'])" \
      "$SRC/firmware-targets.json" "$PROFILE")
IDF_DIR=${ESPSDR_IDF_PATH:-$HOME/.espressif/idf-${PIN:0:8}/esp-idf}
BUILD=${ESPSDR_BUILD_DIR:-$HOME/build/esp32_5g_sniffer/build-$PROFILE}
CMD=${1:-help}
PORT=${2:-/dev/ttyACM0}

activate() {
    [[ -f $IDF_DIR/export.sh ]] || { echo "ESP-IDF not found at $IDF_DIR (see firmware/README.md)"; exit 1; }
    local have
    have=$(git -C "$IDF_DIR" rev-parse HEAD)
    if [[ $have != "$PIN" ]]; then
        echo "ESP-IDF at $IDF_DIR is $have, the firmware pins $PIN"; exit 1
    fi
    # shellcheck disable=SC1091
    . "$IDF_DIR/export.sh" > /dev/null
}

case $CMD in
    env)
        echo "source   $SRC ($(git -C "$SRC" rev-parse --abbrev-ref HEAD) @ $(git -C "$SRC" rev-parse --short HEAD))"
        echo "esp-idf  $IDF_DIR (pinned $PIN)"
        echo "build    $BUILD"
        echo "port     $PORT"
        ;;
    build)
        activate
        mkdir -p "$BUILD"
        # same options as tools/build_firmware.py (the upstream CI), without its export step
        IDF_COMPONENT_MANAGER=0 idf.py -C "$SRC" -B "$BUILD" -DIDF_TARGET=$PROFILE \
            -DSDKCONFIG="$BUILD/sdkconfig" -DSDKCONFIG_DEFAULTS="$SRC/sdkconfig.defaults.$PROFILE" \
            -DESP_SDR_STREAMING=OFF -DRING_PROBE=OFF -DSAMPLE_RATE_PROBE=OFF -DFILTER_REGISTER_PROBE=OFF \
            -DS3_RF_PROBE=OFF -DC5_TUNE_PROBE=OFF -DS2_RF_PROBE=OFF -DC5_GAP_PROBE=OFF build
        echo "built: $BUILD"
        ;;
    flash)
        activate
        IDF_COMPONENT_MANAGER=0 idf.py -C "$SRC" -B "$BUILD" -p "$PORT" flash
        ;;
    erase)
        activate
        python -m esptool --chip esp32c5 -p "$PORT" erase-flash
        ;;
    flash-prebuilt)
        activate
        (cd "$PREBUILT" && sha256sum --ignore-missing -c SHA256SUMS && \
         python -m esptool --chip esp32c5 -p "$PORT" -b 460800 --before default-reset --after hard-reset \
             write-flash @flash_args)
        ;;
    info)
        "$ROOT/.venv_local/bin/python" -m capture.capture_bursts info --port "$PORT" 2>&1 \
            || echo "(run from anywhere; uses $ROOT/.venv_local)"
        ;;
    *)
        sed -n '2,15p' "$0"
        ;;
esac
