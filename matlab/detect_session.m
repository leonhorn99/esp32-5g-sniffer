function T = detect_session(folder, opt)
%DETECT_SESSION  Run ssb_detect on every burst of a capture session and summarise by cell.
%   T = detect_session('data/20261123_110000_n1_2142')
%   T = detect_session(folder, 'fSSB', 2141.05e6, 'Index', 1:50)
%   T = detect_session('data/<survey>', 'FreqMHz', [2610 2625], 'SubcarrierSpacing', 30, 'Lmax', 8, ...
%                      'CFOCenter', 16e3, 'CFOMax', 15e3, 'Tag', 'scs30')    % other bands
%
%   Writes results/<session>_ssb.{md,mat,csv} (derived results never go into data/).
%   T has one row per PSS peak: burst, fLO, fSSB, nid2, pssStat, sssStat, sssOk, pci, iSSB, cfo, snrDb,
%   mibOk, mibOkComb, sfn, tPss_us, tFrame_us.
arguments
    folder (1,:) char
    opt.fSSB double = []                 % default: GSCN search (needs fs >= about 8 MS/s)
    opt.Index double = []
    opt.Verbose (1,1) logical = true
    opt.FreqMHz double = []              % only bursts at these LOs (survey sessions)
    opt.SubcarrierSpacing (1,1) double = 15
    opt.Lmax (1,1) double = 4            % 8 for unpaired (TDD) bands above 2.4 GHz, e.g. n38, n41, n79
    opt.CFOCenter (1,1) double = 0
    opt.CFOMax (1,1) double = 30e3
    opt.Tag (1,:) char = ''              % results/<session>_ssb_<Tag>.*
end
root = fileparts(fileparts(mfilename('fullpath')));
[~, name] = fileparts(strip(folder, 'right', filesep));
% metadata first, then the I/Q in blocks of 200 bursts (a survey session holds GBs of I/Q)
[b, meta] = esp_load_session(folder, 'Index', opt.Index, 'FreqMHz', opt.FreqMHz, 'MetaOnly', true);
if ~isempty(opt.Tag), name = [name '_' opt.Tag]; end
pos = 1:numel(b);
if ~isempty(opt.Index), pos = opt.Index(:).'; end
rows = {};
sssAll = {};                         % 336 SSS statistics per row (empty if the SSB was cut)
tStart = tic;
for c0 = 1:200:numel(b)
    blk = esp_load_session(folder, 'Index', pos(c0:min(end, c0+199)), 'FreqMHz', opt.FreqMHz);
    for kb = 1:numel(blk)
        k = c0 + kb - 1;
        t0 = blk(kb).t0_us;
        if isempty(t0), t0 = []; else, t0 = double(t0); end
        d = ssb_detect(blk(kb).iq, blk(kb).fs, 'fLO', blk(kb).fLO, 'fSSB', opt.fSSB, 't0_us', t0, ...
            'SubcarrierSpacing', opt.SubcarrierSpacing, 'Lmax', opt.Lmax, 'CFOCenter', opt.CFOCenter, 'CFOMax', opt.CFOMax);
        for j = 1:numel(d)
            sfn = NaN;
            if d(j).mibOk, sfn = d(j).mib.sfn; end
            rows(end+1, :) = {blk(kb).idx, blk(kb).fLO/1e6, d(j).fSSB/1e6, d(j).nid2, d(j).pssStat, d(j).sssStat, ...
                d(j).sssOk, d(j).pci, d(j).iSSB, d(j).cfo, d(j).snrDb, d(j).mibOk, d(j).mibOkComb, sfn, ...
                d(j).tPss_us, d(j).tFrame_us}; %#ok<AGROW>
            sssAll{end+1, 1} = d(j).sssAll; %#ok<AGROW>
        end
        if opt.Verbose && mod(k, 25) == 0
            fprintf('  %d/%d bursts, %d PSS peaks, %.0f s\n', k, numel(b), size(rows, 1), toc(tStart));
        end
    end
end
vars = {'burst', 'fLO_MHz', 'fSSB_MHz', 'nid2', 'pssStat', 'sssStat', 'sssOk', 'pci', 'iSSB', 'cfo_Hz', 'snr_dB', ...
    'mibOk', 'mibOkComb', 'sfn', 'tPss_us', 'tFrame_us'};
if isempty(rows)
    T = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
else
    T = cell2table(rows, 'VariableNames', vars);
end

% Cells whose SSS is too weak in single bursts: sum the 336 SSS statistics of all PSS-only
% peaks with the same raster position and N_ID2 (each about Exp(1) under noise, so the sum of
% N is Gamma(N,1)); threshold for P_FA 1e-3 over the 336 hypotheses.
acc = struct('fSSB_MHz', {}, 'nid2', {}, 'n', {}, 'pci', {}, 'stat', {}, 'threshold', {}, ...
    'ratio', {}, 'cfo_Hz', {}, 'tSpread_us', {});
weak = find(~T.sssOk & ~cellfun(@isempty, sssAll));
if ~isempty(weak)
    [g, gf, gn] = findgroups(T.fSSB_MHz(weak), T.nid2(weak));
    for i = 1:max(g)
        idx = weak(g == i);
        if numel(idx) < 3, continue; end
        sm = sum([sssAll{idx}], 2);
        [sv, o] = sort(sm, 'descend');
        thr = gammaincinv(1e-3/336, numel(idx), 'upper');
        tp = T.tPss_us(idx);
        acc(end+1) = struct('fSSB_MHz', gf(i), 'nid2', gn(i), 'n', numel(idx), ...
            'pci', 3*(o(1)-1) + gn(i), 'stat', sv(1), 'threshold', thr, 'ratio', sv(1)/sv(2), ...
            'cfo_Hz', median(T.cfo_Hz(idx)), 'tSpread_us', std(mod(tp - tp(1) + 2500, 5000) - 2500)); %#ok<AGROW>
    end
end

outDir = fullfile(root, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
save(fullfile(outDir, [name '_ssb.mat']), 'T', 'meta', 'acc');
writetable(T, fullfile(outDir, [name '_ssb.csv']));
fid = fopen(fullfile(outDir, [name '_ssb.md']), 'w');
c = onCleanup(@() fclose(fid));
fprintf(fid, '# SSB detection: %s\n\n', name);
fprintf(fid, 'Measured capture, processed by `matlab/detect_session.m` on %s.\n\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
st = meta.settings;
los = unique([b.fLO])/1e6;
fprintf(fid, '- Bursts: %d (CRC ok), LO %s MHz, %g MS/s, %d bit, gain %s\n', numel(b), ...
    mat2str(los), st.rate_msps, st.bits, mat2str(st.gain));
fprintf(fid, '- Search: SSB subcarrier spacing %g kHz, L_max %d, CFO %+.0f +- %.0f kHz, %s\n', ...
    opt.SubcarrierSpacing, opt.Lmax, opt.CFOCenter/1e3, opt.CFOMax/1e3, ...
    ternary(isempty(opt.fSSB), 'all GSCN positions within the capture', sprintf('fSSB %s MHz', mat2str(opt.fSSB/1e6))));
if isfield(meta.device, 'fwopt'), fprintf(fid, '- Firmware: %s\n', meta.device.fwopt); end
rmsv = [b.rms]; clipv = [b.clip_frac];
if ~isempty(rmsv)
    fprintf(fid, '- Burst RMS %.1f counts (median), clipped samples %.3f %% (max)\n', median(rmsv), 100*max(clipv));
end
fprintf(fid, '- PSS peaks above threshold: %d; passed the SSS test: %d; MIB decoded: %d\n\n', ...
    height(T), sum(T.sssOk), sum(T.mibOk | T.mibOkComb));
ok = T(T.sssOk, :);
if isempty(ok)
    fprintf(fid, 'No SSB passed the SSS test.\n');
else
    fprintf(fid, '| SSB centre (MHz) | PCI | SSBs found | SSB indices | MIB decoded | SFNs | CFO median (Hz) | SNR median (dB) |\n');
    fprintf(fid, '|---|---|---|---|---|---|---|---|\n');
    [g, f, p] = findgroups(ok.fSSB_MHz, ok.pci);
    for i = 1:max(g)
        s = ok(g == i, :);
        sfns = unique(s.sfn(~isnan(s.sfn))).';
        fprintf(fid, '| %.2f | %d | %d | %s | %d | %s | %+.0f | %.1f |\n', f(i), p(i), height(s), ...
            mat2str(unique(s.iSSB).'), sum(s.mibOk | s.mibOkComb), mat2str(sfns(1:min(end, 8))), ...
            median(s.cfo_Hz), median(s.snr_dB, 'omitnan'));
    end
end
if ~isempty(acc)
    fprintf(fid, '\n## SSS summed over bursts (cells too weak for single-burst SSS)\n\n');
    fprintf(fid, 'PSS-only peaks with the same raster position and N_ID2; the sum of N SSS statistics is Gamma(N,1) under noise, threshold for P_FA 1e-3 over 336 hypotheses. "timing spread" is the std of the PSS times modulo 5 ms (a real cell gives a small value; noise gives about 1400 us).\n\n');
    fprintf(fid, '| SSB centre (MHz) | N_ID2 | peaks | best PCI | summed stat | threshold | best / 2nd | CFO median (Hz) | timing spread (us) |\n|---|---|---|---|---|---|---|---|---|\n');
    for i = 1:numel(acc)
        a = acc(i);
        pciTxt = sprintf('%d', a.pci);
        if a.stat < a.threshold, pciTxt = sprintf('(%d, below threshold)', a.pci); end
        fprintf(fid, '| %.2f | %d | %d | %s | %.1f | %.1f | %.2f | %+.0f | %.0f |\n', a.fSSB_MHz, a.nid2, a.n, ...
            pciTxt, a.stat, a.threshold, a.ratio, a.cfo_Hz, a.tSpread_us);
    end
end
if opt.Verbose
    for i = 1:numel(acc)
        fprintf('  summed SSS: %.2f MHz N_ID2 %d, %d peaks -> PCI %d (stat %.1f, threshold %.1f)\n', ...
            acc(i).fSSB_MHz, acc(i).nid2, acc(i).n, acc(i).pci, acc(i).stat, acc(i).threshold);
    end
    fprintf('%s: %d bursts, %d PSS peaks, %d SSS ok, %d MIB -> results/%s_ssb.md\n', name, numel(b), ...
        height(T), sum(T.sssOk), sum(T.mibOk | T.mibOkComb), name);
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
