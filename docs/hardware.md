# The ESP32-C5 as a 5G receiver

ESP-SDR (ESPARGOS team) reads raw I/Q from the Wi-Fi/Bluetooth receiver of ESP32 chips through an undocumented debug path; ESPsoup adds long, timestamped bursts and more for the ESP32-C5. This page collects what matters for receiving 5G NR with it. Facts from the firmware projects were checked in November 2026; the firmware is young, so re-check before relying on details.

## Tuning

| | |
|---|---|
| Exact tuning (`FREQ`, whole MHz) | **2130–2730 MHz** and 4790–5745 MHz |
| Tuning path above 3 GHz | ESPsoup uses the Wi-Fi-driver path there by default (`TUNEM 4`, auto): exact up to 5990 MHz. The upstream path would be about 20 dB weaker above 5.75 GHz |
| 1 kHz tuning steps | `FREQK` (ESPsoup) is only approximate (ESPsoup measured +100 kHz → −88 kHz, about ±30 kHz), which would corrupt CFO and crystal measurements. Our tools tune whole MHz and shift digitally, which is exact |
| Below 2130 MHz | the C5 has no 5/6 LO divider: lower `FREQ` values end up near 2.2 GHz |

ESP-SDR's over-the-air CW test of other chips, for comparison: the ESP32, C2, C3, C6, S2 and S3 reach down to 1.85–1.95 GHz (all of band n1), none reached 3.05 GHz. **No ESP32 reaches n78 (3.3–3.8 GHz).**

**5G bands in reach of the C5:**

| Band | Downlink | SSB | Note |
|---|---|---|---|
| n1 | 2110–2170 MHz (FDD) | 15 kHz, case A, L_max 4 | only 2130 MHz and up; an SSB below can still be reached 5–8 MHz below the LO (2125.45 MHz with LO 2131 MHz) |
| n7 | 2620–2690 MHz (FDD) | 15 kHz, case A | NR found here (PCI 317, docs/results.md) |
| n38, n41 | 2496–2690 MHz (TDD) | 15 or 30 kHz, L_max 8 | none found in our survey |
| n40 | 2300–2400 MHz (TDD) | 30 kHz | none found |
| n46, n79 (part) | 5150–5925, 4790–5000 MHz | 30 kHz | none found (not deployed in Europe) |

## Capture modes

- **Bursts:** `CAP16` (8+8 bit) or `CAP20` (10+10 bit) at 80, 40, 20, 10, 8 or 4 MS/s. ESPsoup's `BIGCAP` gives up to 32,704 samples (1.64 ms at 20 MS/s), and `CAPT16/20` adds a timestamp of the first sample (µs since boot).
- **Rates for NR:**
  - 20 MS/s is the safe choice over the air.
  - At 4 and 8 MS/s the samples are simply dropped (no anti-alias filter), so neighbouring carriers fold in.
  - At 40 MS/s the analog filter rolls off beyond about ±8 MHz.
- **Throughput:** USB full speed carries about 0.7–0.9 MB/s, which is 10–13 bursts per second. A continuous 20 MS/s stream is not possible (ESPsoup's `STREAM` stops at 500 kS/s), so SSBDET sends only the SSBs.
- **Gain:** use manual gain (`GAIN MANUAL <index>`); the Wi-Fi AGC settles on packet preambles and clips continuous signals. Gain indices are not monotonic in dB. Auto gain in our tools picks the highest index without clipping (72–80 here), not the best SNR.
- **Data format:** the firmware sends (I, −Q); the host negates Q. Python (`capture/esp_sdr.py`) and MATLAB (`matlab/esp_unpack_iq.m`) decode identically and are tested against the same fixture.

## The crystal

The board's 48 MHz crystal is specified at ±10 ppm, which is up to ±27 kHz at 2.65 GHz, nearly two subcarriers. In practice it is the largest unknown. On one board we measured:

| When | Crystal error | How measured |
|---|---|---|
| first day, cold to warm | −5.1 to −7.3 ppm | CFO of the n1 cell, and SFN timing drift |
| warming up while capturing | −0.2 ppm per minute | live view, 78 s |
| first day, evening, after two spontaneous resets | +2.9 ppm | two cells in n1 and n7 agree within 0.2 ppm |
| second day, right after plugging in → after 8 min of capturing | −3.1 → −6.5 ppm | CFO and SFN timing of the n7 cell |

So the crystal moves by several ppm with temperature, and has jumped across resets. Consequences:
- **Every run must find the CFO first.** The live view scans ±30 kHz with full bursts when started without `--cfo`, and scans again if the board finds no SSB for two runs.
- **The CFO is a measurement of the crystal.** The frequency offset of a cell, and the drift of its SFN-anchored frame timing, measure the crystal against a GNSS-disciplined base station. The two agree to 30–130 Hz, so the chip can calibrate its own clock.

## USB (Linux)

- **Device:** native USB-Serial/JTAG, `/dev/ttyACM0`, ID `303a:1001`; also stable under `/dev/serial/by-id/usb-Espressif_USB_JTAG_serial_debug_unit_<MAC>-if00`.
- **Plugging in:** the firmware often does not answer until the chip is reset once. `capture/ssb_live.py` does this through RTS (`EspSdr.hard_reset()`, nothing is written).
- **Spontaneous resets:** one board did this on one laptop port and not on another; when the board resets, it re-enumerates on USB. The live view reconnects by itself. Possible causes: power, cable, or ModemManager. Firmware version 3 answers `BOOT?` with the reset reason (e.g. `brownout`), and the live view and monitor log it at every reconnect.
- **ModemManager:** tell it to ignore Espressif devices (udev rule in [firmware/README.md](../firmware/README.md)).

## Antenna and levels

- **Antenna:** all results so far used the board's PCB antenna, indoors. Levels are not calibrated: a cable test with a signal generator (`matlab/make_n310_ssb_waveform.m` writes SSB and CW reference files) would give the sensitivity in dBm.
- **What the SSB needs:** MIBs decode down to about −6 dB SNR per resource element. The SSB itself does not stand out in an averaged spectrum: an idle NR cell sends almost only its SSBs (about 1.4 % of the time), so look for cells by correlation, not by spectrum.
