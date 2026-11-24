function test_ssb_chain
%TEST_SSB_CHAIN  End-to-end check of the simulated C5 capture and ssb_detect at high SNR.
%   Run from the project root: matlab -batch "addpath('matlab','tests'); test_ssb_chain"
%   Checks PCI, SSB index, SFN, half-frame bit, CFO (to 150 Hz) and frame timing
%   (to 1 sample) against the generator's ground truth, for the two capture modes of
%   the plan, plus a noise-only burst.
pci = 357;
fSSB = 2140.85e6;            % GSCN 5351 (Magenta's block, 2135-2150 MHz)
% with the 38.211 symbol phase of a real gNB at fSSB: without its compensation in
% ssb_detect the CFO would be off by -fSSB mod 7.0 kHz (= +2.3 kHz here)
[tx, fsTx, tr] = ssb_test_waveform('PCI', pci, 'NumFrames', 4, 'Period', 20, 'CarrierFrequency', fSSB);

% --- mode 1: 20 MS/s, CAP20, 1.64 ms, LO 2138 MHz (SSB 2.85 MHz off DC), full GSCN search
fLO = 2138e6; fs = 20e6; ppm = 8;
t0 = 20e-3 + 0.05e-3;        % burst starts 50 us into frame 2
[rx, m] = c5_rx_model(tx, fsTx, tr.nRE, 'SampleRate', fs, 'N', 32704, ...
    'StartSample', round(t0*fs) + 1, 'FrequencyOffset', fSSB - fLO, ...
    'CarrierFrequency', fSSB, 'PPM', ppm, 'SNRdB', 10, 'Bits', 10, 'Seed', 1, ...
    'DCCounts', 6, 'IQGainDB', 0.3, 'IQPhaseDeg', 2);
tic; [det, info] = ssb_detect(rx, fs, 'fLO', fLO); tSearch = toc;
fprintf('mode 1: %d candidates searched in %.2f s, %d detections\n', numel(info.candidates), tSearch, numel(det));
show(det);
ok = det([det.full] & [det.sssOk]);
inBurst = tr.ssb(tr.ssb.tUseful > m.tStart & tr.ssb.tUseful + 4*(2192/30.72e6) < m.tStart + 32704/fs, :);
assert(numel(ok) == height(inBurst), 'expected %d full SSBs, got %d', height(inBurst), numel(ok));
check(ok, inBurst, m, fSSB, fs);
assert(all([ok.mibOk]), 'MIB not decoded at 10 dB');

% --- mode 2: 4 MS/s, CAP16, 8.18 ms, SSB centred on the LO (cable-test setting)
fs = 4e6; fLO = fSSB;
[rx, m] = c5_rx_model(tx, fsTx, tr.nRE, 'SampleRate', fs, 'N', 32704, ...
    'StartSample', round(19.5e-3*fs) + 1, 'FrequencyOffset', 0, ...
    'CarrierFrequency', fSSB, 'PPM', -6, 'SNRdB', 10, 'Bits', 8, 'Seed', 2, 'DCCounts', 3);
det = ssb_detect(rx, fs, 'fLO', fLO, 'fSSB', fSSB);
fprintf('mode 2: %d detections\n', numel(det));
show(det);
ok = det([det.full] & [det.sssOk]);
inBurst = tr.ssb(tr.ssb.tUseful > m.tStart & tr.ssb.tUseful + 4*(2192/30.72e6) < m.tStart + 32704/fs, :);
assert(numel(ok) == height(inBurst) && height(inBurst) == 4);
check(ok, inBurst, m, fSSB, fs);

% --- noise only, full GSCN search at 20 MS/s
rxn = round(16*(randn(32704,1) + 1i*randn(32704,1))/sqrt(2));
detn = ssb_detect(rxn, 20e6, 'fLO', 2138e6);
assert(isempty(detn) || ~any([detn.sssOk]), 'false SSB detection on noise');
fprintf('noise only: %d PSS peaks above threshold, none passed the SSS test\n', numel(detn));
fprintf('test_ssb_chain: OK\n');
end

function check(ok, inBurst, m, fSSB, fs)
for k = 1:numel(ok)
    d = ok(k);
    tTrue = inBurst.tUseful(k) - m.tStart;            % PSS useful start, from burst start
    assert(d.pci == 357, 'PCI');
    assert(d.iSSB == inBurst.iSSB(k), 'SSB index');
    assert(abs(d.fSSB - fSSB) < 1, 'raster position');
    assert(abs(d.cfo - m.cfoTrue) < 150, 'CFO %.0f Hz vs true %.0f Hz', d.cfo, m.cfoTrue);
    assert(abs(d.tPss - tTrue) <= 1/fs + 1/7.68e6, 'PSS timing %.3g s off', d.tPss - tTrue);
    if d.mibOk
        assert(d.mib.sfn == inBurst.frame(k), 'SFN %d vs %d', d.mib.sfn, inBurst.frame(k));
        assert(d.mib.hrf == inBurst.hrf(k), 'half frame');
        tFrameTrue = inBurst.frame(k)*10e-3 - m.tStart;
        assert(abs(d.tFrame - tFrameTrue) < 2/7.68e6, 'frame timing %.3g s off', d.tFrame - tFrameTrue);
    end
end
end

function show(det)
for d = det(:).'
    fprintf('  f %.2f MHz  nid2 %d  pss %.0f  cfo %+7.0f Hz  t %.4f ms  full %d  sss %.0f  pci %d  iSSB %d  snr %.1f dB  mib %d comb %d', ...
        d.fSSB/1e6, d.nid2, d.pssStat, d.cfo, 1e3*d.tPss, d.full, d.sssStat, d.pci, d.iSSB, d.snrDb, d.mibOk, d.mibOkComb);
    if d.mibOk, fprintf('  SFN %d hrf %d', d.mib.sfn, d.mib.hrf); end
    fprintf('\n');
end
end
