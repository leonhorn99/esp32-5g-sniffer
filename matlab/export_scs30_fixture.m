function export_scs30_fixture()
%EXPORT_SCS30_FIXTURE  30 kHz SSB burst (case C, L_max 8) with ground truth, for tests/test_nr_ssb.py.
%   matlab -batch "addpath('matlab'); export_scs30_fixture"
%   Writes tests/fixtures/ssb_caseC_30k.bin (int16 I,Q interleaved, 15.36 MS/s, the first 2.5 ms
%   of frame 0: all 8 SSBs of half frame 0) and tests/fixtures/ssb_caseC_30k.json (sample rate,
%   PCI, SSB centre frequency, the SSBs with index, half frame, SFN and PSS useful start).
%   The waveform carries the TS 38.211 symbol phase for the RF frequency fSSB, as a gNB sends it.
root = fileparts(fileparts(mfilename('fullpath')));
pci = 501; fSSB = 2595.03e6;          % an n41-like raster position; CarrierFrequency only sets the phase
[wave, fs, tr] = ssb_test_waveform('PCI', pci, 'SubcarrierSpacing', 30, 'BlockPattern', 'Case C', ...
    'TransmittedBlocks', ones(1, 8), 'NumFrames', 1, 'Period', 20, 'SampleRate', 15.36e6, 'CarrierFrequency', fSSB);
n = round(2.5e-3*fs);
x = wave(1:n);
x = x/max(abs([real(x); imag(x)]))*2^14;
v = int16([real(x).'; imag(x).']);
f = fopen(fullfile(root, 'tests', 'fixtures', 'ssb_caseC_30k.bin'), 'w');
fwrite(f, v(:), 'int16');
fclose(f);
t = tr.ssb(tr.ssb.tUseful < 2.5e-3, :);
s = struct('fs', fs, 'pci', pci, 'f_ssb', fSSB, 'scs', 30e3, 'pattern', 'C', 'lmax', 8, ...
    'ssb', struct('issb', num2cell(t.iSSB), 'hrf', num2cell(t.hrf), 'sfn', num2cell(t.frame), ...
                  't_useful_s', num2cell(t.tUseful)));
f = fopen(fullfile(root, 'tests', 'fixtures', 'ssb_caseC_30k.json'), 'w');
fprintf(f, '%s', jsonencode(s));
fclose(f);
fprintf('export_scs30_fixture: %d samples at %.2f MS/s, %d SSBs (PCI %d)\n', n, fs/1e6, height(t), pci);
end
