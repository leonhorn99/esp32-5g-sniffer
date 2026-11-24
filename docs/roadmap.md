# Roadmap: feature ideas

Roughly ordered by value for effort within each group. *(estimate)* marks our guesses.

## Firmware (branch `ssb-sniffer`)

1. **Sub-microsecond timestamps.** Take `t0` from the CPU cycle counter or the systimer's raw ticks instead of `esp_timer` µs. The SFN timing would drop below today's 0.3 µs RMS floor (*estimate:* 20–50 ns). This is the key to time sync and TDOA uses.
2. **On-chip cell search.** All three N_ID2 and several CFO hypotheses (±30 kHz) in acquisition, perhaps over a GSCN raster. The board would find cells without the host's scan, and would recover by itself after the crystal jumps.
3. **SSS and DM-RS on the chip.** Every frame then carries the PCI and the SSB index. Cheap (*estimate:* under 1 ms per SSB).
4. **PBCH/MIB decoding on the chip** (a polar list decoder; *estimate:* 5–20 ms and about 30 KB). This enables a **"cell monitor" mode**: a few bytes per SSB (PCI, SFN, timing, CFO, SNR) instead of 13 KB of I/Q. It works over any link, every SSB of several cells.
5. **Several cells per run.** Round robin over SSB positions, N_ID2 values and frequencies (e.g. n1 and n7). Sectors of one site share the timing.
6. **30 kHz SSBs** (cases B/C) and L_max 8 for TDD bands.
7. **Decimated cuts.** Send the SSB at 4 MS/s (the chip already computes it) instead of raw 20 MS/s: 5× less data, so every SSB of several cells fits over USB. Raw cuts stay optional.
8. **Fine tuning with `FREQK`** (1 kHz steps, ESPsoup), so the LO sits right at the SSB; and **`TUNEM`** above 5.75 GHz.
9. **Robustness.**
   - The firmware should answer after being plugged in without a reset.
   - The reset reason should be reported at start.
   - A watchdog should run around SSBDET.
10. **Wireless sensor node.** Send results over the chip's own Wi-Fi (MQTT/UDP) between captures, from a battery-powered board. *Open question:* the radio is shared between capture and Wi-Fi.

## Host tools

1. **Cell scanner in the GUI.** A band survey with a live spectrum and SSB search over the GSCN raster, a table of cells (frequency, PCI, SNR, MIB), and one click to track a cell.
2. **SIB1 decoding.** The MIB's CORESET#0/searchSpace#0 say where SIB1 is. Capturing those slots and decoding PDCCH and PDSCH would give the operator (PLMN), the global cell ID and the tracking area. A big step (*estimate:* weeks); MATLAB's 5G Toolbox can serve as the reference.
3. **30 kHz and L_max 8 in the NumPy receiver** (today MATLAB only).
4. **Multi-board support.**
   - Several receivers on one cell, aligned by SFN, show the relative timing of the boards: GNSS-free sync, and TDOA experiments.
   - Several cells on one board: relative gNB frequency and timing.
5. **Calibrated levels (RSRP in dBm)** after a cable test with a signal generator, and **coverage maps** with GNSS tags from a phone.
6. **Headless monitor.** A CLI that logs decoded cells to CSV/SQLite and raises alerts for new or unexpected PCIs, frame-timing jumps or frequency outliers. A building block for fake-base-station detection.
7. **Export:** SigMF metadata for captures; a streaming JSON/WebSocket feed of decoded SSBs for other tools.
8. **Crystal tools:** a view of the crystal error and its warm-up; use the cell as a frequency reference to discipline the board's clock.
9. **Channel time series** (Doppler, spectrogram of the per-SSB channel) for sensing experiments.
10. **LTE cell detection** (PSS/SSS/CRS): band 7 and others are full of LTE cells, useful as references and for comparison.
11. **Packaging:** `pip install`, Windows/macOS serial ports, a Jupyter walkthrough of one SSB.

## Measurements worth doing

- **Sensitivity in dBm:** the cable test with a USRP and attenuators (files from `matlab/make_n310_ssb_waveform.m`).
- **Base-station frequency accuracy.** The same board on several cells at the same time cancels the board's crystal. A first hint: two sites differed by about 0.2 ppm in one short run (*estimate*, 3 weak SSBs). 3GPP allows ±0.05 ppm for wide-area base stations.
- **Frame alignment between sites and operators,** and day/night changes, e.g. capacity carriers switched off at night.
- **The crystal's behaviour across resets and temperature,** with a reference thermometer.
