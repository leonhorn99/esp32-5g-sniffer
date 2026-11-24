function out = make_n310_ssb_waveform(opt)
%MAKE_N310_SSB_WAVEFORM  Cable-test files for the N310: SSB-only waveform + CW reference.
%   out = make_n310_ssb_waveform()                       % defaults below
%   out = make_n310_ssb_waveform('fSSB', 2141.05e6, 'PCI', 357, 'OutDir', 'waveforms')
%
%   Writes, for UHD's tx_samples_from_file (--type short = interleaved int16 I,Q):
%     <name>.sc16       SSB burst (case A, 15 kHz, n1-like), NumFrames*10 ms, loops seamlessly
%     <name>_cw.sc16    constant-envelope tone at exactly fSSB, same digital scale as the SSB
%     <name>.json       settings, the N310 command lines and the power relation:
%                       SS-RSRP at the ESP port [dBm] = P_cw [dBm] + rsrpRelCwDb - (extra attenuation)
%                       where P_cw is the CW file's power measured at the same N310 gain and
%                       through the same cable (spectrum analyser or power meter)
%
%   The N310 is tuned to fSSB + LOOffset; the SSB sits at -LOOffset in baseband, so the
%   N310's own LO leakage does not land on the SSB centre. Transmitting in n1 over the air
%   is not allowed: cable plus attenuators only.
arguments
    opt.fSSB (1,1) double = 2141.05e6        % GSCN 5353: LO 2141 MHz puts it 50 kHz off DC
    opt.PCI (1,1) double = 357
    opt.TransmittedBlocks (1,:) double = [1 1 1 1]
    opt.Period (1,1) double = 20
    opt.NumFrames (1,1) double = 8           % file length 80 ms; SFN runs 0..7 and repeats
    opt.SampleRate (1,1) double = 7.68e6     % N310 master clock 122.88 MHz / 16
    opt.LOOffset (1,1) double = -1e6
    opt.Peak (1,1) double = 0.7              % peak |x| relative to int16 full scale
    opt.OutDir (1,:) char = ''
    opt.Addr (1,:) char = '192.168.10.2'     % from the OAI N310 configs (sdr_addrs)
end
root = fileparts(fileparts(mfilename('fullpath')));
if isempty(opt.OutDir), opt.OutDir = fullfile(root, 'waveforms'); end
if ~exist(opt.OutDir, 'dir'), mkdir(opt.OutDir); end

% with the TS 38.211 5.4 symbol phase for the RF frequency of the SSB, as a real gNB sends it
% (the N310 itself upconverts with a continuous-phase LO)
[ssb, fs, tr] = ssb_test_waveform('PCI', opt.PCI, 'TransmittedBlocks', opt.TransmittedBlocks, ...
    'Period', opt.Period, 'NumFrames', opt.NumFrames, 'SampleRate', opt.SampleRate, ...
    'CarrierFrequency', opt.fSSB);
n = (0:numel(ssb)-1).';
fBase = -opt.LOOffset;                        % SSB position in baseband
assert(abs(fBase) + 1.8e6 < fs/2 - 0.2e6, 'SSB does not fit: reduce |LOOffset|');
assert(mod(fBase*numel(ssb)/fs, 1) == 0, 'baseband offset must be periodic over the file');
x = ssb.*exp(1i*2*pi*fBase/fs*n);
s = opt.Peak/max(abs([real(x); imag(x)]));
x = x*s;
cw = opt.Peak*exp(1i*2*pi*fBase/fs*n);       % same RF frequency as the SSB centre

% power relations (per sample, int16 full scale = 1)
Pre = sum(abs(x).^2)*(2048/2192)/tr.nRE/512;  % power of one SSB RE (useful part), Nfft 512
Pcw = opt.Peak^2;
rsrpRelCw = 10*log10(Pre/Pcw);
Pavg = mean(abs(x).^2);

[~, g] = ssb_gscn(opt.fSSB - 1, opt.fSSB + 1);
name = sprintf('ssb_gscn%d_pci%d_p%d_%dblk', g, opt.PCI, opt.Period, sum(opt.TransmittedBlocks));
write_sc16(fullfile(opt.OutDir, [name '.sc16']), x);
write_sc16(fullfile(opt.OutDir, [name '_cw.sc16']), cw);

fTx = opt.fSSB + opt.LOOffset;
uhd = @(file) sprintf(['tx_samples_from_file --args "addr=%s,master_clock_rate=122.88e6" ' ...
    '--file %s --type short --rate %.0f --freq %.0f --gain <GAIN> --subdev "A:0" --ant "TX/RX" --repeat'], ...
    opt.Addr, file, fs, fTx);
out = struct('name', name, 'fSSB_Hz', opt.fSSB, 'gscn', g, 'pci', opt.PCI, ...
    'nid2', mod(opt.PCI, 3), 'blockPattern', tr.pattern, 'transmittedBlocks', opt.TransmittedBlocks, ...
    'period_ms', opt.Period, 'fileLength_ms', 10*opt.NumFrames, 'sampleRate', fs, ...
    'n310Freq_Hz', fTx, 'basebandOffset_Hz', fBase, 'peak', opt.Peak, ...
    'rsrpRelCwDb', rsrpRelCw, 'avgPowerRelCwDb', 10*log10(Pavg/Pcw), ...
    'ssbBandPowerRelRsrpDb', 10*log10(240), ...
    'cmdSSB', uhd([name '.sc16']), 'cmdCW', uhd([name '_cw.sc16']), ...
    'espLoSuggestion_MHz', struct('rate4MSps', round(opt.fSSB/1e6), 'rate20MSps', round(opt.fSSB/1e6) - 3), ...
    'note', ['SS-RSRP at the ESP port = P_cw(dBm, measured) + rsrpRelCwDb - extra attenuation (dB). ' ...
             'UHD flags: check tx_samples_from_file --help of the installed UHD version.']);
fid = fopen(fullfile(opt.OutDir, [name '.json']), 'w');
fwrite(fid, jsonencode(out, 'PrettyPrint', true));
fclose(fid);
fprintf('%s: %d samples at %.2f MS/s (%.0f ms), RSRP = P_cw %+.1f dB\n', name, numel(x), fs/1e6, ...
    1e3*numel(x)/fs, rsrpRelCw);
fprintf('  %s\n  %s\n', out.cmdSSB, out.cmdCW);
end

function write_sc16(file, x)
v = zeros(2, numel(x), 'int16');
v(1, :) = int16(round(32767*real(x)));
v(2, :) = int16(round(32767*imag(x)));
fid = fopen(file, 'w', 'ieee-le');
fwrite(fid, v(:), 'int16');
fclose(fid);
end
