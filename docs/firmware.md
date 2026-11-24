# Firmware

**Where the firmware lives:**
- The board runs our fork of ESPsoup, branch `ssb-sniffer`, in the submodule `firmware/espsoup-firmware`.
- Building, flashing, the pinned ESP-IDF and USB quirks: [firmware/README.md](../firmware/README.md).
- The full SSBDET protocol and algorithm: [firmware/espsoup-firmware/docs/ssbdet.md](../firmware/espsoup-firmware/docs/ssbdet.md).

What our branch adds to ESPsoup 0.7.1 is one command, `SSBDET`: on-chip detection of the SSBs of one NR cell. The board sends only the raw I/Q around each SSB, with its timestamp.

## Using SSBDET from the host

Command line (stores the cuts as a session in `data/`):

```bash
.venv_local/bin/python -m capture.capture_bursts ssb-stream --port /dev/ttyACM0 \
    --freq 2131 --ssb 2125.45 --nid2 0 --cfo 12000 --gain 76 --duration 30 --label mycell
```

Live view (decodes and shows every SSB; `--record` stores them):

```bash
.venv_local/bin/python -m capture.ssb_live --port /dev/ttyACM0 --freq 2131 --ssb 2125.45 [--nid2 0 --cfo 12000]
```

Library (`capture/esp_sdr.py`):

```python
from capture.esp_sdr import EspSdr, ssb_stream
sdr = EspSdr("/dev/ttyACM0"); sdr.connect(); sdr.set_freq_mhz(2131); sdr.set_gain(76)
for kind, item in ssb_stream(sdr, rate_code=2, offset_hz=-5_550_000, nid2=0, cfo_hz=12_000, duration=10):
    if kind == "frame":       # SsbFrame: .iq(), .t_first_us, .pss_idx - .start, .cfo_hz, .stat, .mode
        print(item.seq, item.mode, item.stat, item.cfo_hz)
    else:                     # "start", "end", "idle" (alive, no SSB), "text" (board log lines)
        print(kind, item)
```

**What the parameters mean:**
- `--nid2` is the PCI mod 3.
- `--cfo` is the expected frequency offset, i.e. the board's crystal error × frequency. It must be within about ±2 kHz; the board then tracks it.
- The live view finds both by itself when they are left out.

## Performance and limits

| | |
|---|---|
| SSB cuts per second (20 MS/s, 8 bit, 20 ms cell) | about 25, i.e. every SSB; limited by USB full speed |
| cut length | 0.33 ms: the SSB plus 20 µs on either side, 13 KB at 8 bit |
| tracking misses | 7 in 743 over 30 s |
| on-chip time per SSB | decimation about 5 ms, PSS search 1–2 ms (tracking), send about 18 ms |
| timestamp | 1 µs resolution (`esp_timer`), which sets the 0.3 µs RMS timing floor |

**Limits:**
- One cell per run, 15 kHz SSB spacing only.
- The cell's N_ID2 and a rough CFO must be known; the host scans for them.

Ideas for the firmware are in [roadmap.md](roadmap.md).

## Testing the firmware DSP on a PC

`ssb_dsp.h` has no ESP-IDF dependency, so `firmware/espsoup-firmware/tests/ssb/host_ssb.c` compiles it on the PC. `tests/test_ssb_onchip_dsp.py` runs it on recorded bursts (the sample `samples/n1_2131_20M_30bursts`, or the full captures) and compares with the MATLAB detector:
- 18/18 SSBs found;
- positions within 0.17 µs;
- no false detection.

Run it after every change to `ssb_dsp.h` and before flashing:

```bash
.venv_local/bin/python -m unittest tests.test_ssb_onchip_dsp -v
```
