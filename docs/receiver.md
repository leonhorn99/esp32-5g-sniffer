# The NR SSB receiver

Two implementations of the same chain:
- `matlab/ssb_detect.m`: the reference. It is batch only and uses the 5G Toolbox. It covers the GSCN search, 15/30 kHz and L_max 4/8.
- `capture/nr_ssb.py`: NumPy, used live: 15 and 30 kHz SSBs, L_max 4 and 8 (cases A, B, C), the PSS search over a whole raster, RSRP.

The firmware does the first steps on the chip: the PSS (all versions) and the SSS → PCI (version 3).

## The chain

1. **Mix and search.**
   - Shift the candidate SSB centre (SS_REF, subcarrier 120 of the SSB) to 0 Hz.
   - Correlate with the three PSS over a grid of CFO hypotheses, after a ±1.1 MHz low-pass.
   - The statistic is normalised to a unit-mean exponential under white noise for one trial. The thresholds are 18 for a burst search and 14 on the chip (10.5 in its ±40 µs tracking window).
2. **CFO inside the PSS symbol:** the phase between the correlations of its two halves.
3. **Fine timing:**
   - MATLAB: the correlation peak at 7.68 MS/s.
   - Python: at the raw rate with parabolic interpolation (sub-sample).
4. **OFDM demodulation of the four SSB symbols.**
   - The FFT window starts 9 samples (at 7.68 MS/s) into the cyclic prefix; that ramp is undone.
   - **The TS 38.211 §5.4 symbol phase is removed** (see below).
5. **SSS → PCI.**
   - The SSS is correlated against all 336 N_ID1, coherently with the PSS channel (smoothed over 9 subcarriers). The threshold is 13: 336 · e⁻¹³ ≈ 10⁻³ false PCIs per PSS peak.
   - For cells too weak for one burst, `detect_session` sums the 336 statistics of all PSS peaks at one raster position (Gamma(N, 1) under noise).
6. **CFO refinement from the PSS → SSS phase** (two symbols apart). It is unambiguous only within ±3.5 kHz (period 7.0 kHz). If the MIB fails the CRC, the two candidates one period away are tried, and one is kept only if its MIB passes.
7. **PBCH DM-RS → SSB index and half frame (L_max 4).**
   - MATLAB: coherent against the PSS/SSS channel.
   - Python: coherent within 12-subcarrier chunks across the three symbols.
8. **Channel estimate** from all pilots (PSS, SSS, DM-RS), smoothed over ±6 subcarriers. The noise variance comes from the DM-RS residuals.
9. **PBCH:**
   - LLRs from the equalised QPSK symbols;
   - descrambling with c(i + v·864), v = i_SSB mod 4;
   - rate recovery (the 864 bits repeat the 512-bit code);
   - **CRC-aided polar SC-list decoding** (list 8).
10. **BCH:** CRC24C, payload descrambling, de-interleaving (TS 38.212 §7.1) → the 24-bit MIB, the 4 SFN LSBs, the half-frame bit and the k_SSB MSB.
11. **Timing:** the PSS time minus the SSB's symbol position in the half frame, minus 5 ms × half frame, gives the frame start. Frame start minus SFN × 10 ms is constant on the cell's clock.

## Searching a whole raster quickly

`nr_ssb.pss_search_raster` searches all GSCN positions within a burst (about 30 within ±6.1 MHz at 20 MS/s), all three N_ID2 and the CFO hypotheses at once:
- **One FFT of the burst with 250 Hz bins.** GSCN offsets (multiples of 50 kHz) and CFO steps (5 kHz) are then whole bins, so each hypothesis is just a ±1.1 MHz slice of the same spectrum.
- **One small inverse FFT per hypothesis**, which correlates at 5 MS/s. At 2.5 MS/s the 0.4 µs lag steps let the peak fall between lags and lost up to 2.7 dB; at 5 MS/s the loss against the full-rate search is 0.13 dB median, 0.59 dB at most.
- **Speed:** about 0.3 s per 1.64 ms burst for 30 positions × 7 CFOs × 3 N_ID2, decoding included.
- **On the 23 Nov captures:** 18 of 18 SSBs of PCI 165, no false cell; on the n7 survey bursts it finds PCI 317 like MATLAB.

The cell scanner (`capture/cell_scan.py`) and the live view's `--ssb auto` use it.

## 30 kHz SSBs and L_max 8

`receive(..., scs=30e3, lmax=8, pattern="C")`:
- **Symbol timing:** symbols are half as long, and the long CP is 144 + 16·2^µ samples at 2048·SCS, on the first symbol of every 0.5 ms.
- **SSB index:** for L_max 8 the DM-RS carries all 3 bits of the SSB index; the half frame comes from the MIB.
- **PBCH scrambling:** uses i_SSB mod 8.
- **Validation:** checked on a case C burst from `nrWaveformGenerator` with the RF symbol phase (`matlab/export_scs30_fixture.m`, 8 SSBs). All 8 give the right PCI and SSB index, the MIB, the CFO within 2 Hz and the frame start to 0.05 µs.

## RSRP

`receive()` reports `rsrp_dbfs`: the mean channel power on the SSS subcarriers (SS-RSRP is defined on the SSS resource elements, TS 38.215), in dB relative to the 10-bit ADC full scale per resource element. Payload units are converted for CAP16 and for decimated cuts. It is not calibrated: the gain index adds an unknown, frequency-dependent gain, and a cable test with a known level is needed for dBm.

## Pitfalls found on real signals

### The TS 38.211 §5.4 symbol phase (found 23 Nov 2026)

TS 38.211 §5.4 defines the upconversion per OFDM symbol: every subcarrier's phase restarts at the useful start of every symbol, at its absolute RF frequency.

**What a receiver sees:** tuned to the SSB, it sees a phase step of 2π·f_SSB·T_sym from one symbol to the next. A phase-based CFO estimate across symbols reads that step as an extra frequency offset of −f_SSB mod 7.0 kHz.

**On real signals:**
- **+3.93 kHz at 2125.45 MHz.** We measured +3.89 and +3.99 kHz above the crystal error, and first took it for a synthesizer offset.
- **The DM-RS search lost coherence across symbols**, so SSB indices came out wrong on weak SSBs.

**Why the tests missed it:** synthetic signals from `nrWaveformGenerator` carry this phase only if its `CarrierFrequency` is set. With the default 0 they have none, so the simulation could not show the problem.

**Checked in MATLAB's own OFDM model:** `nrOFDMModulate` with `CarrierFrequency`, the SSB placed at three positions in a 52-RB carrier, at three carrier frequencies.
- In all 9 cases the apparent CFO is exactly −f_SSB mod 7007.3 Hz.
- It does **not** depend on the carrier centre f0.
- Multiplying symbol l by exp(j2π f_SSB l T_sym) removes it completely.

**Fix:** both receivers now do that, with the absolute SSB frequency, and the synthetic waveforms carry the phase (`ssb_test_waveform(..., 'CarrierFrequency', fSSB)`). On the recorded sessions:
- the CFO equals the crystal error from the frame timing, to 30–130 Hz;
- the SNR rose by 1.3–1.5 dB;
- the MIBs of a weak sector went from 4/5 to 5/5, and on one session from 597 to 724.

### The 7.0 kHz ambiguity

The PSS → SSS refinement spans two symbols, so its period is 1/(2 T_sym) = 7.007 kHz. The PSS-halves estimate before it can be off by more than half a period at low SNR. The CRC rule above resolves it. Choosing the candidate by the DM-RS metric alone made the AWGN simulation slightly worse.

### The CFO search window

The C5's crystal moves by several ppm (warm-up, resets; [hardware.md](hardware.md)). A search window centred on 0 Hz can miss the cell, while a strong cell is still found with a 9 kHz error and looks fine. So:
- centre the search on the crystal error measured on a known cell (`CFOCenter` = −ppm · f);
- or scan ±30 kHz (the live view does).

### Timing floor

The firmware's `t0_us` has 1 µs resolution. Uniform rounding alone gives 1/√12 = 0.29 µs RMS, and that is what the SFN-anchored timing reaches (0.29–0.30 µs RMS over 30–78 s). A phase-slope timing estimate across subcarriers was noisier (0.51 µs) and was replaced by the interpolated correlation peak. Firmware version 3 stamps in 62.5 ns steps; how far that lowers the floor on the chip is not yet measured.

## Validation

- **Bit-exact against the 5G Toolbox** (`tests/test_nr_ssb.py`, vectors from `matlab/export_pbch_tables.m`):
  - PSS, SSS, PBCH DM-RS, the PBCH and DM-RS indices, `nrBCH` and `nrPBCH`, for six cells;
  - decoding back to the MIB;
  - the polar information positions and rate matching are derived from `nrPolarEncode`/`nrRateMatchPolar`, not typed in.
- **Noise:**
  - 20/20 MIBs at −2 dB per coded bit;
  - random LLRs never pass the CRC.
- **Against the MATLAB detector on 150 recorded cuts:**
  - 148 MIBs against 147, with the same PCI, SFN and SSB index;
  - CFO within +4 Hz (median), PSS time within 0.036 µs (spread).
- **MATLAB end to end on simulated C5 captures** (`tests/test_ssb_chain.m`):
  - PCI, SSB index, SFN, half frame, CFO to 150 Hz and timing to one sample, at 20 MS/s and 4 MS/s;
  - no false SSB on noise.
- **Physical consistency on real data:**
  - the SFN advances with the capture timestamps;
  - the frame timing is constant to 0.3 µs;
  - the CFO matches the crystal error measured independently from the timing;
  - two base stations agree on the board's clock within 0.2 ppm;
  - re-encoding a decoded MIB matches 82 % of the received coded bits, against 45 % for a wrong one.
