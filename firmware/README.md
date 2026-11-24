# Firmware for the ESP32-C5

The board runs **ESPsoup**, the ESP32-C5 fork of the ESPARGOS team's [ESP-SDR](https://github.com/ESPARGOS/esp-sdr) (GPL-3.0-or-later), with our additions on the branch **`ssb-sniffer`**. It is receive-only firmware:
- `CAPT16`/`CAPT20` bursts of up to 32,704 samples with timestamps, manual gain and DC calibration (from ESPsoup);
- `SSBDET`, our on-chip SSB detection, which sends only the raw I/Q around each detected SSB.

| Path | What |
|---|---|
| `espsoup-firmware/` | git submodule: our fork of [pit711/espsoup-firmware](https://github.com/pit711/espsoup-firmware), branch `ssb-sniffer` (based on ESPsoup tag v0.7.1). Its own submodule `components/esp-dsp` (v1.8.2) |
| `fw.sh` | build, flash, flash an official prebuilt image, ask the board for its info |
| `prebuilt/` | (not in git) put an official ESPsoup release image here for `fw.sh flash-prebuilt` |

SSBDET itself (command, frame format, algorithm, timing) is documented in [espsoup-firmware/docs/ssbdet.md](espsoup-firmware/docs/ssbdet.md); the host side is in [../docs/firmware.md](../docs/firmware.md).

## Boards

Tested on one ESP32-C5 development board (ESP32-C5HR8, chip revision v1.0, 48 MHz crystal, 16 MB flash, 8 MB PSRAM) over its native USB (USB-Serial/JTAG, `/dev/ttyACM0` on Linux). The firmware's 2 MB flash setting works on larger flash.

The factory demo on such boards does not use the USB serial port. If you want to keep it, back up the whole flash first (`esptool --chip esp32c5 -p /dev/ttyACM0 read-flash 0 ALL factory.bin`). `fw.sh erase` before the first flash also removes the demo's stored PHY calibration.

## ESP-IDF: a pinned version

The firmware reads the receiver through undocumented functions of Espressif's binary PHY libraries, and these change between ESP-IDF versions. Its `firmware-targets.json` therefore pins the C5 build to ESP-IDF commit `25fe69f946311abdaf9ad56591f25fedbc20ac98` (a master commit of 2026-07-11, reporting v6.2.0), and its CI refuses any other.

- Install that commit next to any other ESP-IDF, e.g. as a shallow clone in `~/.espressif/idf-25fe69f9/esp-idf`, then run its `install.sh esp32c5`.
- `fw.sh` expects it there; set `ESPSDR_IDF_PATH` (and `ESPSDR_BUILD_DIR` for the build folder) to use other locations.
- `fw.sh build` checks the commit before building.

Our build of ESPsoup's tag v0.7.1 with this commit is byte-identical in code to the official 0.7.1 release image. It differs only in the app descriptor (version string, build time, IDF name, ELF hash) and the final checksum.

## Build and flash

```bash
firmware/fw.sh build                  # -> ~/build/esp32_5g_sniffer/build-esp32c5 (outside the source tree)
firmware/fw.sh flash /dev/ttyACM0     # our build
firmware/fw.sh flash-prebuilt         # an official ESPsoup image from firmware/prebuilt/ instead
firmware/fw.sh info                   # INFO / CAPS / LIMITS / VERSION / FWOPT from the board
```

The build uses the same options as the upstream CI (`tools/build_firmware.py`): target esp32c5, `sdkconfig.defaults.esp32c5`, streaming and all diagnostic probes off, component manager off. After flashing, `fw.sh info` should answer `C5SDR 6 burst 16380`, list `BIGCAP CAPT ... SSBDET` in CAPS, and report FWOPT `0.7.1-dev+ssb2`.

### In VS Code with the ESP-IDF extension

The extension treats the opened workspace folder as the ESP-IDF project and has no setting for a subfolder. With the repo root open, *Set Espressif Device Target* fails with "CMakeLists.txt not found". Therefore:

1. Open `ESP32_5G_Sniffer.code-workspace` (*File → Open Workspace from File…*). It holds the repo root and `firmware (ESPsoup C5)`.
2. In the ESP-IDF status bar, pick the folder `firmware (ESPsoup C5)` (*ESP-IDF: Pick a Workspace Folder*).
3. Select the pinned ESP-IDF for that folder, with build directory `~/build/esp32_5g_sniffer/build-esp32c5` and `sdkconfig.defaults.esp32c5`. That is the same build as `fw.sh`.

Without the extension: *Terminal → Run Task* offers `firmware: build`, `firmware: flash`, `firmware: flash prebuilt 0.7.1` and `firmware: info` (they call `fw.sh`).

## USB quirks (Linux)

- **Flashing** over the C5's native USB resets the chip into the bootloader by itself. If not: hold BOOT, tap RESET, release BOOT.
- **After plugging the board in**, the firmware often does not answer until the chip is reset once.
  - `capture/ssb_live.py` does this by itself (an RTS reset through the USB-Serial/JTAG port; nothing is written).
  - By hand: `esptool --chip esp32c5 -p /dev/ttyACM0 --before default-reset --after hard-reset read-mac`.
- **ModemManager** probes new serial ports and can reset or confuse the board. Tell it to ignore Espressif devices with `/etc/udev/rules.d/99-espressif-no-modemmanager.rules` containing `ATTRS{idVendor}=="303a", ENV{ID_MM_DEVICE_IGNORE}="1"`, then `sudo udevadm control --reload-rules`.

## Provenance

- ESPsoup's release manifest for 0.7.1 names source commit `2dcf034c`, which is not in the public history. The public tag v0.7.1 is `085c752e` (same day), and our branch starts from it.
- Upstream ESP-SDR has fixes newer than ESPsoup, e.g. "Fix capture clock stalls on C5 and C61" (2026-11-16). They are candidates to cherry-pick onto `ssb-sniffer` after checking them on hardware.
- Licence: GPL-3.0-or-later, like ESP-SDR and ESPsoup.
