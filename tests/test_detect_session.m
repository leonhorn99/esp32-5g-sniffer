function test_detect_session
%TEST_DETECT_SESSION  detect_session on a simulated capture session (session format of capture/session.py).
%   Run from the project root: matlab -batch "addpath('matlab','tests'); test_detect_session"
%   Six CAP16 bursts at 20 MS/s, LO 2139 MHz, SSB of PCI 101 at 2141.05 MHz (GSCN 5353),
%   SNR 0 dB per RE: bursts 1-3 hold SSBs, bursts 4-6 hold none. Full GSCN search.
root = fileparts(fileparts(mfilename('fullpath')));
fs = 20e6; fLO = 2139e6; fSSB = 2141.05e6; N = 32704;
[tx, fsTx, tr] = ssb_test_waveform('PCI', 101, 'NumFrames', 4, 'Period', 20, 'CarrierFrequency', fSSB);
starts = [20.05e-3, 20.10e-3, 0.08e-3, 5e-3, 10e-3, 14e-3];   % first three contain SSBs
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder, 's'));
fid = fopen(fullfile(folder, 'bursts.bin'), 'w');
jl = fopen(fullfile(folder, 'bursts.jsonl'), 'w');
offset = 0;
for k = 1:numel(starts)
    rx = c5_rx_model(tx, fsTx, tr.nRE, 'SampleRate', fs, 'N', N, 'StartSample', round(starts(k)*fs) + 1, ...
        'FrequencyOffset', fSSB - fLO, 'CarrierFrequency', fSSB, 'PPM', 5, 'SNRdB', 0, 'Bits', 8, 'Seed', k);
    v = int8([real(rx).'; -imag(rx).']);                  % device order: I, -Q
    fwrite(fid, v(:), 'int8');
    rec = struct('idx', k-1, 'offset', offset, 'nbytes', 2*N, 'fmt', 'CAP16', 'n', N, 'rate_code', 2, ...
        'fs', fs, 'freq_mhz', fLO/1e6, 'gain', 50, 'bandwidth_mhz', [], 'crc_ok', true, 'crc32', '00000000', ...
        'capture_us', 821, 't0_us', 1e6 + 1e6*starts(k), 'seq', k-1, 'lag', 7, 'host_time', 0, ...
        'temp_c', [], 'reply', '', 'rms', 16, 'clip_frac', 0);
    fprintf(jl, '%s\n', jsonencode(rec));
    offset = offset + 2*N;
end
fclose(fid); fclose(jl);
meta = struct('format', 1, 'device', struct('info', 'C5SDR 6 burst 16380', 'fwopt', 'simulated'), ...
    'settings', struct('freq_mhz', fLO/1e6, 'rate_msps', 20, 'bits', 8, 'gain', 50));
fid = fopen(fullfile(folder, 'session.json'), 'w'); fwrite(fid, jsonencode(meta)); fclose(fid);

T = detect_session(folder, 'Verbose', false);
[~, name] = fileparts(folder);
outMd = fullfile(root, 'results', [name '_ssb.md']);
outMat = fullfile(root, 'results', [name '_ssb.mat']);
c2 = onCleanup(@() delete(outMd, outMat));
ok = T(T.sssOk, :);
assert(~isempty(ok) && all(ok.pci == 101), 'PCI 101 not found or wrong PCI');
assert(all(abs(ok.fSSB_MHz - 2141.05) < 1e-6), 'wrong raster position');
assert(all(ismember(ok.burst, 0:2)), 'detection in a burst without SSB');
assert(any(ok.mibOk | ok.mibOkComb), 'no MIB decoded at 0 dB');
assert(all(ok.sfn(~isnan(ok.sfn)) == [2*ones(nnz(~isnan(ok.sfn(ok.burst < 2))), 1); zeros(nnz(~isnan(ok.sfn(ok.burst == 2))), 1)]), 'SFN');
assert(all(abs(ok.cfo_Hz - (-5e-6*fSSB)) < 300), 'CFO');
txt = fileread(outMd);
assert(contains(txt, '| 2141.05 | 101 |'), 'summary table row missing');
fprintf('test_detect_session: OK (%d SSBs found in bursts %s, %d with MIB)\n', height(ok), ...
    mat2str(unique(ok.burst).'), sum(ok.mibOk | ok.mibOkComb));
end
