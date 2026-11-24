function [wave, fs, truth] = ssb_test_waveform(opt)
%SSB_TEST_WAVEFORM  SSB-only NR downlink waveform with SS_REF at 0 Hz, plus ground truth.
%   [wave, fs, truth] = ssb_test_waveform('PCI', 357, 'NumFrames', 2, ...)
%   The carrier is exactly the 20 RBs of the SSB (NCRBSSB = 0, KSSB = 0), so the
%   SSB centre (subcarrier 120) sits at 0 Hz. The MIB comes from nrWaveformGenerator.
%   CarrierFrequency: the RF frequency the SSB centre will be received at (Hz). It adds the
%   TS 38.211 5.4 symbol phase a real gNB applies (0 = none; ssb_detect expects it).
%
%   truth.ssb  table: frame (SFN), hrf, iSSB, symbol, tUseful (s, start of the PSS
%              useful part from the waveform start)
%   truth.nRE  number of SSB resource elements in the whole waveform
%   truth.reMagnitudes  unique RE magnitudes in the grid (all equal is assumed by the
%              SNR definition in c5_rx_model)
arguments
    opt.PCI (1,1) double = 357
    opt.SubcarrierSpacing (1,1) double = 15
    opt.BlockPattern (1,:) char = ''
    opt.TransmittedBlocks (1,:) double = [1 1 1 1]
    opt.Period (1,1) double = 20
    opt.NumFrames (1,1) double = 2
    opt.SampleRate (1,1) double = 7.68e6
    opt.CarrierFrequency (1,1) double = 0
end
scs = opt.SubcarrierSpacing;
pattern = opt.BlockPattern;
if isempty(pattern)
    pattern = 'Case A';
    if scs == 30, pattern = 'Case C'; end
end

cfg = nrDLCarrierConfig;
cfg.FrequencyRange = 'FR1';
cfg.ChannelBandwidth = 5*(scs/15);   % 5 MHz (15 kHz) or 10 MHz (30 kHz): both allow 20 RBs
cfg.NCellID = opt.PCI;
cfg.NumSubframes = 10*opt.NumFrames;
cfg.InitialNSubframe = 0;
cfg.SampleRate = opt.SampleRate;
cfg.CarrierFrequency = opt.CarrierFrequency;   % carrier centre = SS_REF here
sc = nrSCSCarrierConfig('SubcarrierSpacing', scs, 'NSizeGrid', 20, 'NStartGrid', 0);
cfg.SCSCarriers = {sc};
bwp = nrWavegenBWPConfig('SubcarrierSpacing', scs, 'NSizeBWP', 20, 'NStartBWP', 0);
cfg.BandwidthParts = {bwp};
ssb = nrWavegenSSBurstConfig;
ssb.BlockPattern = pattern;
ssb.TransmittedBlocks = opt.TransmittedBlocks;
ssb.Period = opt.Period;
ssb.NCRBSSB = 0;
ssb.KSSB = 0;
ssb.SubcarrierSpacingCommon = scs;
cfg.SSBurst = ssb;
% SSB only: switch off the default PDSCH, PDCCH and CSI-RS
for k = 1:numel(cfg.PDSCH), cfg.PDSCH{k}.Enable = false; end
for k = 1:numel(cfg.PDCCH), cfg.PDCCH{k}.Enable = false; end
for k = 1:numel(cfg.CSIRS), cfg.CSIRS{k}.Enable = false; end

[wave, info] = nrWaveformGenerator(cfg);
fs = info.ResourceGrids(1).Info.SampleRate;
% The SSB is not part of info.ResourceGrids: demodulate the waveform to get its REs.
carrier = nrCarrierConfig('SubcarrierSpacing', scs, 'NSizeGrid', 20, 'NStartGrid', 0);
grid = nrOFDMDemodulate(carrier, wave, 'SampleRate', fs);
mag = abs(grid(:));
nz = mag(mag > 1e-3*max(mag));

% ground truth: SSB first symbols in a half frame
L = numel(opt.TransmittedBlocks);
symIdx = ssb_first_symbols(pattern, L);
rows = [];
halfFramePeriod = opt.Period/5;             % in half frames
for hfAbs = 0:(2*opt.NumFrames - 1)
    if mod(hfAbs, halfFramePeriod) ~= 0, continue; end
    for i = 0:L-1
        if ~opt.TransmittedBlocks(i+1), continue; end
        t = hfAbs*5e-3 + ssb_symbol_time(symIdx(i+1), scs);
        rows = [rows; floor(hfAbs/2), mod(hfAbs, 2), i, symIdx(i+1), t]; %#ok<AGROW>
    end
end
truth.ssb = array2table(rows, 'VariableNames', {'frame', 'hrf', 'iSSB', 'symbol', 'tUseful'});
truth.nRE = numel(nz);
truth.reMagnitudes = uniquetol(nz, 1e-6);
truth.pci = opt.PCI;
truth.scs = scs;
truth.pattern = pattern;
truth.Lmax = L;
truth.fs = fs;
truth.cfg = cfg;
end

function s = ssb_first_symbols(pattern, L)
switch pattern
    case {'Case A', 'Case C'}
        n = 0:(L/2 - 1);
        s = reshape([2; 8] + 14*n, 1, []);
    case 'Case B'
        n = 0:(L/4 - 1);
        s = reshape([4; 8; 16; 20] + 28*n, 1, []);
    otherwise
        error('ssb_test_waveform:pattern', 'pattern %s not supported', pattern);
end
end
