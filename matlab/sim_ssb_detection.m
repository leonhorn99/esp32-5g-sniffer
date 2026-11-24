function R = sim_ssb_detection(opt)
%SIM_SSB_DETECTION  Monte Carlo of SSB detection and MIB decoding for the C5 capture modes.
%   R = sim_ssb_detection()                    % full run, writes results/sim_ssb_detection.*
%   R = sim_ssb_detection('Trials', 5, 'SNRdB', [-6 0], 'NoiseBursts', 3, 'Tag', 'smoke')
%
%   SIMULATION, not a measurement. Channel: AWGN, one path, static. Receiver model:
%   c5_rx_model (CFO from a +-10 ppm crystal, DC, IQ imbalance, 8/10-bit quantisation;
%   no aliasing of the 4/8 MS/s modes, no AGC, no interference). SNR is per resource
%   element (SS-SINR style). The detector searches only the true raster position in
%   the SNR sweeps; the noise-only bursts of mode A use the full GSCN search.
%
%   Modes (from the plan):
%     A  20 MS/s, CAP20 (10 bit), 32,704 samples = 1.64 ms, SSB 3.05 MHz above the LO
%     B   4 MS/s, CAP16 (8 bit),  32,704 samples = 8.18 ms, SSB 50 kHz from the LO
arguments
    opt.Trials (1,1) double = 100
    opt.SNRdB (1,:) double = -14:2:4
    opt.NoiseBursts (1,1) double = 100
    opt.Seed (1,1) double = 2026
    opt.Tag (1,:) char = ''
end
rng(opt.Seed);
root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
name = 'sim_ssb_detection';
if ~isempty(opt.Tag), name = [name '_' opt.Tag]; end

fSSB = 2141.05e6;                                   % GSCN 5353
pcis = [0 1 2 101 357 512 763 1007];
W = struct([]);
for k = 1:numel(pcis)
    [tx, fsTx, tr] = ssb_test_waveform('PCI', pcis(k), 'NumFrames', 8, 'Period', 20, 'CarrierFrequency', fSSB);
    W(k).tx = tx; W(k).tr = tr;
    W(k).Ere = sum(abs(tx).^2)/fsTx*(2048/2192)/tr.nRE;
end

modes = struct('name', {'A: 20 MS/s, 10 bit, 1.64 ms', 'B: 4 MS/s, 8 bit, 8.18 ms'}, ...
    'fs', {20e6, 4e6}, 'bits', {10, 8}, 'fOff', {3.05e6, 0.05e6}, ...
    'lead', {[0.02 0.12]*1e-3, [0.05 6.3]*1e-3});
N = 32704;
fsS = 7.68e6;
ssbDur = 3*548/fsS + 512/fsS;                       % PSS useful start to the end of symbol 3
nS = numel(opt.SNRdB);
for im = 1:numel(modes)
    md = modes(im);
    acc = zeros(nS, 9);  % ssbTotal, pssFound, pciOk, mibOk, bursts, burstPci, burstMib, burstMibOrComb, falsePci
    cfoErr = cell(nS, 1); snrErr = cell(nS, 1);
    tStartClock = tic;
    for is = 1:nS
        for t = 1:opt.Trials
            w = W(randi(numel(W)));
            frame = 2*randi(3);                     % SS bursts at frames 0,2,4,6; use 2..6
            ssbT = w.tr.ssb(w.tr.ssb.frame == frame, :);
            tStart = ssbT.tUseful(1) - (md.lead(1) + rand*diff(md.lead));
            margin = 64/fsTx;
            i0 = floor((tStart - margin)*fsTx) + 1;
            i1 = ceil((tStart + N/md.fs + margin)*fsTx) + 1;
            seg = w.tx(i0:i1);
            s0 = round((tStart - (i0-1)/fsTx)*md.fs) + 1;
            tStart = (i0-1)/fsTx + (s0-1)/md.fs;
            [rx, m] = c5_rx_model(seg, fsTx, w.tr.nRE, 'EnergyPerRE', w.Ere, 'SampleRate', md.fs, ...
                'N', N, 'StartSample', s0, 'FrequencyOffset', md.fOff, 'CarrierFrequency', fSSB, ...
                'PPM', 20*rand - 10, 'SNRdB', opt.SNRdB(is), 'Bits', md.bits, ...
                'DCCounts', 4*rand, 'IQGainDB', rand - 0.5, 'IQPhaseDeg', 6*rand - 3);
            [det, info] = ssb_detect(rx, md.fs, 'fLO', fSSB - md.fOff, 'fSSB', fSSB);
            tTrue = ssbT.tUseful - tStart;
            full = tTrue - 9/fsS >= 0 & tTrue + ssbDur <= N/md.fs;
            matched = false(size(det));
            burstPci = false; burstMib = false;
            for k = find(full).'
                j = find(abs([det.tPss] - tTrue(k)) < 2/fsS, 1);
                acc(is, 1) = acc(is, 1) + 1;
                if isempty(j), continue; end
                matched(j) = true;
                d = det(j);
                acc(is, 2) = acc(is, 2) + 1;
                pciOk = d.sssOk && d.pci == w.tr.pci;
                mibOk = pciOk && d.mibOk && d.mib.sfn == frame && d.mib.hrf == ssbT.hrf(k);
                acc(is, 3) = acc(is, 3) + pciOk;
                acc(is, 4) = acc(is, 4) + mibOk;
                burstPci = burstPci || pciOk;
                burstMib = burstMib || mibOk;
                if pciOk
                    cfoErr{is}(end+1) = d.cfo - m.cfoTrue;
                    snrErr{is}(end+1) = d.snrDb - opt.SNRdB(is);
                end
            end
            combOk = any(arrayfun(@(c) c.mibOk && c.pci == w.tr.pci && c.mib.sfn == frame, info.combined));
            if ~isempty(det)
                bad = [det.sssOk] & (~matched | [det.pci] ~= w.tr.pci);
                acc(is, 9) = acc(is, 9) + sum(bad);
            end
            acc(is, 5:8) = acc(is, 5:8) + [1, burstPci, burstMib, burstMib || combOk];
        end
        fprintf('%s  SNR %+5.1f dB: P(PCI | SSB) %.2f  P(MIB | SSB) %.2f  P(MIB or combined | burst) %.2f  (%.0f s)\n', ...
            md.name, opt.SNRdB(is), acc(is,3)/acc(is,1), acc(is,4)/acc(is,1), acc(is,8)/acc(is,5), toc(tStartClock));
    end
    R.mode(im).name = md.name;
    R.mode(im).snrDb = opt.SNRdB;
    R.mode(im).counts = acc;
    R.mode(im).countNames = {'ssbTotal', 'pssFound', 'pciOk', 'mibOk', 'bursts', 'burstPci', 'burstMib', 'burstMibOrComb', 'falsePci'};
    R.mode(im).pPss = acc(:,2)./acc(:,1);
    R.mode(im).pPci = acc(:,3)./acc(:,1);
    R.mode(im).pMib = acc(:,4)./acc(:,1);
    R.mode(im).pBurstMibComb = acc(:,8)./acc(:,5);
    R.mode(im).ssbPerBurst = acc(:,1)./acc(:,5);
    R.mode(im).cfoRmsHz = cellfun(@(v) sqrt(mean(v.^2)), cfoErr);
    R.mode(im).snrEstBiasDb = cellfun(@mean, snrErr);
end

% noise only: false alarms. Mode A with the full GSCN search around a 2138 MHz LO.
fa = zeros(2, 3);   % bursts, PSS peaks, SSS passes
for b = 1:opt.NoiseBursts
    rxn = round(64*(randn(N,1) + 1i*randn(N,1))/sqrt(2));
    det = ssb_detect(rxn, 20e6, 'fLO', 2138e6);
    fa(1,:) = fa(1,:) + [1, numel(det), sum([det.sssOk])];
    rxn = round(16*(randn(N,1) + 1i*randn(N,1))/sqrt(2));
    det = ssb_detect(rxn, 4e6, 'fLO', fSSB - 0.05e6, 'fSSB', fSSB);
    fa(2,:) = fa(2,:) + [1, numel(det), sum([det.sssOk])];
end
[~, nCand] = ssb_gscn(2138e6 - (10e6 - 1.9e6), 2138e6 + (10e6 - 1.9e6));
R.noise = array2table(fa, 'VariableNames', {'bursts', 'pssPeaks', 'sssPasses'}, ...
    'RowNames', {sprintf('A, full GSCN search (%d positions)', numel(nCand)), 'B, one position'});
R.opt = opt;
R.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
save(fullfile(outDir, [name '.mat']), 'R');
write_report(R, fullfile(outDir, [name '.md']));
plot_sim_ssb_detection(R, fullfile(outDir, [name '.png']));
fprintf('wrote %s.{mat,md,png} in %s\n', name, outDir);
end

function write_report(R, file)
fid = fopen(file, 'w');
c = onCleanup(@() fclose(fid));
fprintf(fid, '# SSB detection with the simulated ESP32-C5 capture (simulation)\n\n');
fprintf(fid, 'Written by `matlab/sim_ssb_detection.m` on %s. %d bursts per SNR point and mode.\n', R.created, R.opt.Trials);
fprintf(fid, 'AWGN, static single path, no fading, no interference, no aliasing; CFO from a +-10 ppm crystal, DC, IQ imbalance, quantisation. ');
fprintf(fid, 'SNR is per resource element (SS-SINR style). Not measured.\n\n');
for im = 1:numel(R.mode)
    M = R.mode(im);
    fprintf(fid, '## %s\n\n', M.name);
    fprintf(fid, 'Full SSBs per burst: %.2f on average (they were placed inside the burst).\n\n', mean(M.ssbPerBurst));
    fprintf(fid, '| SNR per RE (dB) | P(PSS found) | P(PCI correct) | P(MIB, one SSB) | P(MIB in burst, incl. combining) | CFO error RMS (Hz) | SNR estimate bias (dB) | false PCIs |\n');
    fprintf(fid, '|---|---|---|---|---|---|---|---|\n');
    for is = 1:numel(M.snrDb)
        fprintf(fid, '| %+.0f | %.2f | %.2f | %.2f | %.2f | %.0f | %+.1f | %d |\n', M.snrDb(is), M.pPss(is), ...
            M.pPci(is), M.pMib(is), M.pBurstMibComb(is), M.cfoRmsHz(is), M.snrEstBiasDb(is), M.counts(is, 9));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '## Noise only (false alarms)\n\n| Case | bursts | PSS peaks above threshold | passed the SSS test |\n|---|---|---|---|\n');
for k = 1:height(R.noise)
    fprintf(fid, '| %s | %d | %d | %d |\n', R.noise.Properties.RowNames{k}, R.noise.bursts(k), ...
        R.noise.pssPeaks(k), R.noise.sssPasses(k));
end
end

