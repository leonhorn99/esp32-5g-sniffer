function res = sib1_offline(folder, opt)
%SIB1_OFFLINE  Decode SIB1 of one NR cell from a recorded session of full bursts (offline).
%   res = sib1_offline('data/20261123_165332_survey_2175_2715', 'FreqMHz', 2670, 'SSBMHz', 2667.75, 'PCI', 317)
%
%   A burst (1.64 ms at 20 MS/s) that holds the SSB never holds the slots of SIB1's control channel
%   (2-5 ms later), so the steps are:
%   1. SSBs: ssb_detect on every burst at the LO; the bursts with this PCI's MIB give its frame start
%      on the board's clock (tFrame_us) and the SFN. A line through them (tFrame_us = a + b x frame)
%      puts every burst on the cell's frame grid.
%   2. From the MIB (TS 38.213 V17.14.0, Tables 13-1 and 13-11, as capture/nr_mib.py): CORESET#0
%      (RBs, symbols, offset below the SSB's resource-block grid, k_SSB) and the Type0-PDCCH slots
%      n0, n0+1 of even or odd frames. The carrier grid here is CORESET#0 itself (NStartGrid 0, the
%      reference point of CORESET#0's PDCCH and of the SIB1 PDSCH's DM-RS).
%   3. Every slot n0 / n0+1 that lies whole in a burst: mix the CORESET#0 centre to 0 Hz (with the
%      SSB's CFO scaled to that frequency), resample to 1024 x 15 kHz, OFDM-demodulate with the
%      TS 38.211 5.4 symbol phase at the absolute frequency, blind-decode the Type0 common search
%      space (AL 4, 8, 16) with the SI-RNTI (nrPDCCHDecode, nrDCIDecode), over timing offsets
%      opt.Timing_us around the predicted slot start.
%   4. DCI format 1_0 with SI-RNTI (TS 38.212 V17.13.0 7.3.1.2.1): FDRA (RIV over N_RB of
%      CORESET#0), TDRA (TS 38.214 V17.16.0 Table 5.1.2.1.1-2, default A), VRB-to-PRB, MCS (Table
%      5.1.3.1-1), RV, SI indicator. PDSCH as TS 38.214 5.1.6.2 has it for DCI 1_0: DM-RS type 1,
%      single symbol, pos2, CDM groups without data 2 (1 for 2-symbol type B), port 1000; VRB
%      bundles of 2. Then nrPDSCHDecode and nrDLSCHDecoder: a passed CRC is a decoded SIB1.
%
%   SI messages other than SIB1 (DCI with SI indicator 1, Type0A search space) differ in two ways, both
%   needing the decoded SIB1 (capture/nr_sib1.py --matlab prints the options):
%   - the PDSCH DM-RS sequence starts at point A (CRB 0), not at CORESET#0 ('OffsetToPointA':
%     CORESET#0 starts at CRB offsetToPointA - Table 13-1 offset);
%   - the time allocation comes from SIB1's pdsch-TimeDomainAllocationList when it has one
%     ('TimeAllocationList'), not from default table A (TS 38.214 Table 5.1.2.1.1-1).
%   Found on PCI 317 on 26 Nov: with both, all 8 copies of its SI message decoded; without, none.
%   5. 'Combine', true (weak cells): (a) every slot whose DCI decoded but whose PDSCH did not is
%      demodulated with its own DCI (allocation, RV) and added to a DL-SCH soft buffer per transport
%      block size until the CRC passes (PCI 165 moves its SIB1 to other RBs each time); (b) if all decoded
%      DCIs have the same allocation (PCI 317), it is also assumed in the slots without a decoded DCI
%      whose PDSCH DM-RS is there (metric >= PresenceThreshold). NR signals the RV only in the DCI (TS 38.321 5.3.1 has no rule for SIB1); the RV of an
%      occasion is taken as ceil(3k/2) mod 4 with k = (SFN/2) mod 4 (the LTE rule, TS 36.321 5.3.1),
%      which every DCI decoded so far followed. A passing CRC confirms it.
%
%   res  struct: anchors, the timing fit, MIB-derived geometry, and .tries (one row per slot tried:
%        burst, SFN, slot, timing offset, DCI found, DCI fields, TBS, CRC, SIB1 bits as hex)
%   Writes results/sib1_<session>_<PCI>.mat and .json (decoded SIB1s as hex, for the ASN.1 step).
arguments
    folder (1,:) char
    opt.FreqMHz (1,1) double
    opt.SSBMHz (1,1) double
    opt.PCI (1,1) double
    opt.CFOCenter (1,1) double = 0
    opt.CFOMax (1,1) double = 30e3
    opt.Timing_us double = [0 -1 1 -2 2 -3 3 -4 4 -5 5]
    opt.Out (1,:) char = ''
    opt.KeepGrids (1,1) logical = false            % keep each tried slot's grid (at its best DM-RS offset) in res.tries
    opt.Combine (1,1) logical = false              % soft-combine the PDSCH of every occasion (see below)
    opt.CombineSlot double = []                    % slot to combine (default: the slot of the first decoded DCI)
    opt.PresenceThreshold (1,1) double = 0.35      % PDSCH DM-RS metric for 'a PDSCH is there' (noise ~0.17)
    opt.OffsetToPointA (1,1) double = NaN          % from the decoded SIB1: needed for other SI messages (SIB2...)
    opt.TimeAllocationList double = []             % SIB1 pdsch-TimeDomainAllocationList: rows [typeA(1)/B(0) S L]
    opt.Slots double = []                          % try these slots of the frames below instead of SIB1's n0, n0+1
    opt.FrameMod (1,1) double = 1                  % ... in frames with mod(SFN, FrameMod) in FrameSet
    opt.FrameSet double = 0
end
t0clock = tic;
[~, name] = fileparts(folder);
bursts = esp_load_session(folder, 'FreqMHz', opt.FreqMHz);
fprintf('%s: %d bursts at LO %g MHz\n', name, numel(bursts), opt.FreqMHz);
fSSB = opt.SSBMHz * 1e6;

% 1. anchors: SSBs of this PCI with a decoded MIB
anc = struct('burst', {}, 'tFrame_us', {}, 'sfn', {}, 'cfo', {}, 'iSSB', {}, 'mib', {}, 'snrDb', {});
for b = 1:numel(bursts)
    x = bursts(b);
    det = ssb_detect(x.iq, x.fs, 'fLO', x.fLO, 'fSSB', fSSB, 'CFOCenter', opt.CFOCenter, 'CFOMax', opt.CFOMax, ...
                     't0_us', double(x.t0_us));
    for d = det(:)'
        if d.mibOk && d.pci == opt.PCI
            anc(end + 1) = struct('burst', b, 'tFrame_us', d.tFrame_us, 'sfn', d.mib.sfn, 'cfo', d.cfo, ...
                                  'iSSB', d.iSSB, 'mib', d.mib, 'snrDb', d.snrDb); %#ok<AGROW>
            break
        end
    end
end
fprintf('anchors (PCI %d with MIB): %d\n', opt.PCI, numel(anc));
if numel(anc) < 2
    error('sib1_offline:anchors', 'need at least two bursts with this PCI''s MIB for the frame clock');
end
% frame numbers: SFN unwrapped with the board-time spacing (10.24 s per SFN cycle)
[~, o] = sort([anc.tFrame_us]);
anc = anc(o);
k = zeros(1, numel(anc));
k(1) = anc(1).sfn;
for j = 2:numel(anc)
    approx = k(1) + (anc(j).tFrame_us - anc(1).tFrame_us) / 1e4;
    k(j) = anc(j).sfn + 1024 * round((approx - anc(j).sfn) / 1024);
end
p = polyfit(k, [anc.tFrame_us], 1);              % tFrame_us = p(2) + p(1) * frame
resid = [anc.tFrame_us] - polyval(p, k);
fprintf('frame clock: %.6f us per frame (board vs cell %+.3f ppm), residual RMS %.3f us, max %.3f us\n', ...
        p(1), (p(1) / 1e4 - 1) * 1e6, rms(resid), max(abs(resid)));
mib = anc(1).mib;
cfoSSB = median([anc.cfo]);
iSSB = mode([anc.iSSB]);

% 2. CORESET#0 and Type0-PDCCH occasions from the MIB
[c0, ss0] = sib1_tables();
if mib.kSSB >= 24
    error('sib1_offline:noSIB1', 'k_SSB %d >= 24: no CORESET for SIB1', mib.kSSB);
end
if mib.scsCommonFlag ~= 0
    error('sib1_offline:scs', 'subCarrierSpacingCommon 30 kHz: Table 13-2 is not included');
end
row = c0(bitshift(mib.pdcchConfigSIB1, -4) + 1, :);   % [pattern nRB nSym offset]
nRB = row(2); nSymC = row(3); offRB = row(4);
s0 = ss0(bitand(mib.pdcchConfigSIB1, 15) + 1, :);      % [O perSlot M firstSym(-1 = alt)]
v = s0(1) + floor(iSSB * s0(3));
n0 = mod(v, 10); parity = mod(floor(v / 10), 2);
evenOdd = {'even', 'odd'};
if s0(4) < 0
    firstSym = (mod(iSSB, 2) == 1) * nSymC;
else
    firstSym = s0(4);
end
scs = 15e3;
fLow = fSSB - 120 * scs - mib.kSSB * scs - offRB * 12 * scs;   % subcarrier 0 of CORESET#0
f0 = fLow + nRB * 6 * scs;                                     % DC of a grid of nRB RBs starting there
fprintf(['MIB: pdcch-ConfigSIB1 %d -> CORESET#0 %d RB x %d sym, offset %d RB, %.3f-%.3f MHz; ' ...
         'SIB1 PDCCH in slots %d,%d of %s frames from symbol %d (SSB index %d); CFO at the SSB %+.0f Hz\n'], ...
        mib.pdcchConfigSIB1, nRB, nSymC, offRB, fLow / 1e6, (fLow + nRB * 12 * scs) / 1e6, n0, n0 + 1, ...
        evenOdd{parity + 1}, firstSym, iSSB, cfoSSB);

carrier = nrCarrierConfig('SubcarrierSpacing', 15, 'NSizeGrid', nRB, 'NStartGrid', 0, 'NCellID', opt.PCI);
cs = nrCORESETConfig('CORESETID', 0, 'FrequencyResources', ones(1, nRB / 6), 'Duration', nSymC, ...
                     'CCEREGMapping', 'interleaved', 'REGBundleSize', 6, 'InterleaverSize', 2, 'ShiftIndex', opt.PCI);
nCCE = nRB * nSymC / 6;
ss = nrSearchSpaceConfig('SearchSpaceID', 0, 'CORESETID', 0, 'SearchSpaceType', 'common', ...
                         'StartSymbolWithinSlot', firstSym, 'SlotPeriodAndOffset', [1 0], 'Duration', 1, ...
                         'NumCandidates', [0 0 min(4, floor(nCCE / 4)) min(2, floor(nCCE / 8)) min(1, floor(nCCE / 16))]);
pdcch = nrPDCCHConfig('CORESET', cs, 'SearchSpace', ss, 'NStartBWP', 0, 'NSizeBWP', nRB, 'RNTI', 0, ...
                      'DMRSScramblingID', opt.PCI);  % common search space: no RNTI hashing (SI-RNTI is checked in
                                                 % nrDCIDecode); DM-RS N_ID = PCI (the toolbox's default is 2)
nFDRA = ceil(log2(nRB * (nRB + 1) / 2));
siOpt = struct('c0crb', opt.OffsetToPointA - offRB, 'tdra', opt.TimeAllocationList);   % for SI messages (see help)
Kdci = nFDRA + 4 + 1 + 5 + 2 + 1 + 15;               % DCI 1_0, SI-RNTI, licensed spectrum

% 3. slots n0, n0+1 of the right frames that lie whole in a burst
fsNew = 15.36e6 * 2 ^ (bursts(1).fs > 2 * 15.36e6);  % 1024 x 15 kHz (2048 for 40 MS/s)
cand = zeros(0, 4);                                   % burst, frame, slot, start sample (0-based, fractional)
for b = 1:numel(bursts)
    x = bursts(b);
    dur_us = numel(x.iq) / x.fs * 1e6;
    fA = (double(x.t0_us) - p(2)) / p(1);
    for F = floor(fA) - 1:ceil(fA + dur_us / 1e4) + 1
        if isempty(opt.Slots)
            if mod(mod(F, 1024), 2) ~= parity
                continue
            end
            slotList = [n0 n0 + 1];
        else
            if ~ismember(mod(mod(F, 1024), opt.FrameMod), opt.FrameSet)
                continue
            end
            slotList = opt.Slots;
        end
        for s = slotList
            ts = p(2) + p(1) * (F + s / 10);
            st = (ts - double(x.t0_us)) * x.fs / 1e6;
            if st - 6 * x.fs / 1e6 >= 0 && st + (1e3 + 6) * x.fs / 1e6 <= numel(x.iq)
                cand(end + 1, :) = [b F s st]; %#ok<AGROW>
            end
        end
    end
end
fprintf('slots %d,%d of %s frames lying whole in a burst: %d\n', n0, n0 + 1, evenOdd{parity + 1}, size(cand, 1));

nfft = round(fsNew / scs);
okFail = {'fail', 'ok'};
tries = struct('burst', {}, 'sfn', {}, 'slot', {}, 'dt_us', {}, 'dci', {}, 'al', {}, 'cand', {}, 'fields', {}, ...
               'tbs', {}, 'crcOk', {}, 'hex', {}, 'evmPdcch', {}, 'snrPdsch', {}, 'dmrsMetric', {}, ...
               'dmrsDt_us', {}, 'symPowDb', {}, 'grid', {}, 'grid0', {});
lastB = 0;
for c = 1:size(cand, 1)
    b = cand(c, 1); F = cand(c, 2); s = cand(c, 3);
    x = bursts(b);
    if b ~= lastB
        n = (0:numel(x.iq) - 1).';
        cfo0 = cfoSSB * f0 / fSSB;
        y = x.iq .* exp(-2j * pi * (f0 - x.fLO + cfo0) * n / x.fs);
        [pp, qq] = rat(fsNew / x.fs);
        y = resample(y, pp, qq);
        lastB = b;
    end
    carrier.NSlot = s;
    carrier.NFrame = mod(F, 1024);
    t = struct('burst', b, 'sfn', mod(F, 1024), 'slot', s, 'dt_us', NaN, 'dci', false, 'al', NaN, 'cand', NaN, ...
               'fields', [], 'tbs', NaN, 'crcOk', false, 'hex', '', 'evmPdcch', NaN, 'snrPdsch', NaN, ...
               'dmrsMetric', NaN, 'dmrsDt_us', NaN, 'symPowDb', [], 'grid', [], 'grid0', []);
    for dt = opt.Timing_us
        st = round((cand(c, 4) + dt * x.fs / 1e6) * fsNew / x.fs);
        seg = y(st + 1:st + round(fsNew * 1e-3));
        grid = nrOFDMDemodulate(carrier, seg, 'SampleRate', fsNew, 'Nfft', nfft, 'CarrierFrequency', f0);
        if opt.Combine && dt == 0
            t.grid0 = grid;                        % the frame clock's own timing (0.1-0.7 us RMS)
        end
        dm = sib1_dmrs_metric(carrier, pdcch, grid, nCCE);
        if ~(dm <= t.dmrsMetric)                   % best over the timing offsets (NaN at first)
            t.dmrsMetric = dm; t.dmrsDt_us = dt;
            pw = mean(abs(grid) .^ 2, 1);
            t.symPowDb = 10 * log10(pw / median(pw));
            if opt.KeepGrids
                t.grid = grid;
            end
        end
        [dcibits, al, ci, evm] = sib1_blind_pdcch(carrier, pdcch, grid, Kdci, opt.PCI);
        if ~isempty(dcibits) && opt.Combine
            t.grid0 = grid;                        % the timing whose PDCCH decoded
        end
        if ~isempty(dcibits)
            t.dt_us = dt; t.dci = true; t.al = al; t.cand = ci; t.evmPdcch = evm;
            t.fields = sib1_dci_fields(dcibits, nFDRA, nRB);
            [t.tbs, t.crcOk, blk, t.snrPdsch] = sib1_pdsch(carrier, grid, t.fields, mib.dmrsTypeAPosition, ...
                                                           opt.PCI, nRB, siOpt);
            if t.crcOk
                t.hex = lower(reshape(dec2hex(bin2dec(reshape(char(blk(:).' + '0'), 8, []).'), 2).', 1, []));
            end
            break
        end
    end
    tries(end + 1) = t; %#ok<AGROW>
    if t.dci
        f = t.fields;
        fprintf(['burst %3d SFN %4d slot %d dt %+g us: DCI AL%d c%d (EVM %.0f%%): RB %d+%d, TDRA %d, ' ...
                 'VRB %d, MCS %d, RV %d, SI %d -> TBS %d, CRC %s\n'], b, t.sfn, s, t.dt_us, t.al, t.cand, ...
                100 * t.evmPdcch, f.rbStart, f.nRB, f.tdra, f.vrbInterleaved, f.mcs, f.rv, f.si, t.tbs, ...
                okFail{t.crcOk + 1});
    else
        fprintf('burst %3d SFN %4d slot %d: no DCI with SI-RNTI; CORESET DM-RS metric %.2f (dt %+g us); symbol power dB %s\n', ...
                b, t.sfn, s, t.dmrsMetric, t.dmrsDt_us, mat2str(round(t.symPowDb, 1)));
    end
end

comb = struct('used', 0, 'decoded', {{}});
if opt.Combine
    comb = sib1_combine(tries, carrier, mib, opt, nRB, okFail, siOpt);
    for k = 1:numel(comb.decoded)
        d = comb.decoded{k};
        tries(end + 1) = struct('burst', NaN, 'sfn', d.sfn, 'slot', d.slot, 'dt_us', 0, 'dci', false, 'al', NaN, ...
            'cand', NaN, 'fields', d.fields, 'tbs', d.tbs, 'crcOk', true, 'hex', d.hex, 'evmPdcch', NaN, ...
            'snrPdsch', NaN, 'dmrsMetric', NaN, 'dmrsDt_us', NaN, 'symPowDb', [], 'grid', [], 'grid0', []); %#ok<AGROW>
    end
    [tries.grid0] = deal([]);                      % the grids are not needed in the result
end

res = struct('session', name, 'pci', opt.PCI, 'fLO', opt.FreqMHz * 1e6, 'fSSB', fSSB, 'anchors', anc, ...
             'clock', p, 'clockResidRms_us', rms(resid), 'mib', mib, 'iSSB', iSSB, 'cfoSSB', cfoSSB, ...
             'coreset0', struct('nRB', nRB, 'nSym', nSymC, 'offsetRB', offRB, 'fLow', fLow, 'f0', f0), ...
             'type0', struct('n0', n0, 'frameParity', parity, 'firstSymbol', firstSym), 'tries', tries, ...
             'combine', rmfield(comb, 'decoded'));
fprintf('SIB1: %d slots tried, %d DCIs, %d CRC ok (%.0f s)\n', numel(tries), sum([tries.dci]), sum([tries.crcOk]), ...
        toc(t0clock));
out = opt.Out;
if isempty(out)
    out = fullfile(fileparts(mfilename('fullpath')), '..', 'results', sprintf('sib1_%s_%d', name, opt.PCI));
end
save([out '.mat'], 'res');
ok = tries([tries.crcOk]);
j = struct('session', name, 'pci', opt.PCI, 'ssb_mhz', opt.SSBMHz, 'lo_mhz', opt.FreqMHz, 'kssb', mib.kSSB, ...
           'pdcch_config_sib1', mib.pdcchConfigSIB1, 'ssb_index', iSSB, ...
           'sib1', struct('sfn', {ok.sfn}, 'slot', {ok.slot}, 'tbs', {ok.tbs}, 'hex', {ok.hex}, ...
                          'si', num2cell(arrayfun(@(t) t.fields.si, ok))));   % si 0: SIB1, 1: other SI message
fid = fopen([out '.json'], 'w');
fprintf(fid, '%s\n', jsonencode(j, 'PrettyPrint', true));
fclose(fid);
end


function [c0, ss0] = sib1_tables()
% TS 38.213 V17.14.0 (ETSI TS 138 213, 2026-09), checked against the ETSI PDF on 24 Nov 2026; the same
% tables as capture/nr_mib.py.
% Table 13-1 ({SS/PBCH block, PDCCH} SCS {15, 15} kHz, min. channel bandwidth 5 or 10 MHz):
% [pattern, N_RB, N_symb, offset in RBs]; index 15 reserved (NaN)
c0 = [1 24 2 0; 1 24 2 2; 1 24 2 4; 1 24 3 0; 1 24 3 2; 1 24 3 4; 1 48 1 12; 1 48 1 16; 1 48 2 12; ...
      1 48 2 16; 1 48 3 12; 1 48 3 16; 1 96 1 38; 1 96 2 38; 1 96 3 38; NaN NaN NaN NaN];
% Table 13-11 (pattern 1, FR1): [O, sets per slot, M, first symbol (-1: 0 for even i, N_symb^CORESET for odd i)]
ss0 = [0 1 1 0; 0 2 1/2 -1; 2 1 1 0; 2 2 1/2 -1; 5 1 1 0; 5 2 1/2 -1; 7 1 1 0; 7 2 1/2 -1; ...
       0 1 2 0; 5 1 2 0; 0 1 1 1; 0 1 1 2; 2 1 1 1; 2 1 1 2; 5 1 1 1; 5 1 1 2];
end


function [dcibits, al, ci, evm] = sib1_blind_pdcch(carrier, pdcch, grid, K, ncellid)
% Type0 common search space: every candidate of AL 4, 8, 16; the first whose CRC matches the SI-RNTI.
dcibits = []; al = NaN; ci = NaN; evm = NaN;
nc = pdcch.SearchSpace.NumCandidates;
for L = [4 8 16]
    for c = 1:nc(log2(L) + 1)
        pdcch.AggregationLevel = L;
        pdcch.AllocatedCandidate = c;
        [ind, dmrs, dmrsInd] = nrPDCCHResources(carrier, pdcch);
        if isempty(ind)
            continue
        end
        [h, nVar] = nrChannelEstimate(carrier, grid, dmrsInd, dmrs);
        [sym, hs] = nrExtractResources(ind, grid, h);
        eq = nrEqualizeMMSE(sym, hs, nVar);
        cw = nrPDCCHDecode(eq, ncellid, 0, nVar);
        [bits, mask] = nrDCIDecode(cw, K, 8, 65535);
        if mask == 0
            dcibits = bits; al = L; ci = c;
            ref = nrPDCCH(nrDCIEncode(bits, 65535, L * 108), ncellid, 0);
            evm = sqrt(mean(abs(eq - ref) .^ 2) / mean(abs(ref) .^ 2));
            return
        end
    end
end
end


function m = sib1_dmrs_metric(carrier, pdcch, grid, nCCE)
% How well the PDCCH DM-RS of the whole CORESET matches, per REG bundle: |sum conj(d) r|^2 /
% (N sum |r|^2), about 1/N for noise and up to 1 for a clean PDCCH; the mean of the best half of the
% bundles (a PDCCH of AL 4 fills half of an 8-CCE CORESET).
L = 2 ^ floor(log2(nCCE));
pdcch.AggregationLevel = min(L, 16);
nc = [0 0 0 0 0];
nc(log2(pdcch.AggregationLevel) + 1) = 1;
pdcch.SearchSpace.NumCandidates = nc;
pdcch.AllocatedCandidate = 1;
[~, d, di] = nrPDCCHResources(carrier, pdcch);
r = grid(di);
nb = numel(d) / 18;                                 % 6 REGs x 3 DM-RS per bundle
c = reshape(conj(d) .* r, 18, nb);
e = reshape(abs(r) .^ 2, 18, nb);
q = sort(abs(sum(c, 1)) .^ 2 ./ (18 * sum(e, 1)), 'descend');
m = mean(q(1:max(1, floor(nb / 2))));
end


function f = sib1_dci_fields(b, nFDRA, nRB)
% DCI format 1_0 with CRC scrambled by SI-RNTI (TS 38.212 7.3.1.2.1), MSB first.
v = @(i, n) bin2dec(char(double(b(i:i + n - 1)).' + '0'));
riv = v(1, nFDRA); p = nFDRA + 1;
f.tdra = v(p, 4); p = p + 4;
f.vrbInterleaved = v(p, 1); p = p + 1;
f.mcs = v(p, 5); p = p + 5;
f.rv = v(p, 2); p = p + 2;
f.si = v(p, 1);
f.riv = riv;
% RIV (TS 38.214 5.1.2.2.2): L - 1 <= floor(N/2): RIV = N(L-1) + start, else N(N-L+1) + (N-1-start)
L = floor(riv / nRB) + 1; st = mod(riv, nRB);
if L + st > nRB
    L = nRB - L + 2; st = nRB - 1 - st;
end
f.rbStart = st; f.nRB = L;
end


function [tbs, ok, blk, snrDb] = sib1_pdsch(carrier, grid, f, dmrsPos, ncellid, nRB, siOpt)
% Decode the PDSCH a DCI 1_0 with SI-RNTI schedules, from one slot.
if nargin < 7
    siOpt = struct('c0crb', NaN, 'tdra', []);
end
[pdsch, R] = sib1_pdsch_config(f, dmrsPos, ncellid, nRB, siOpt);
if f.si == 1 && ~isnan(siOpt.c0crb)               % DM-RS from point A: CORESET#0 at its CRB index
    carrier.NStartGrid = siOpt.c0crb;
    pdsch.NStartBWP = siOpt.c0crb;
end
[cw, snrDb, ~, tbs] = sib1_pdsch_llr(carrier, grid, pdsch, R);
dec = nrDLSCHDecoder('TargetCodeRate', R, 'TransportBlockLength', tbs, 'LDPCDecodingAlgorithm', 'Normalized min-sum', ...
                     'MaximumLDPCIterationCount', 20);
[blk, err] = dec(cw, pdsch.Modulation, 1, f.rv);
ok = ~err;
end


function [pdsch, R] = sib1_pdsch_config(f, dmrsPos, ncellid, nRB, siOpt)
% The PDSCH of a DCI 1_0 with SI-RNTI (TS 38.214 V17.16.0 5.1.2.1.1 Table 5.1.2.1.1-2, 5.1.3.1, 5.1.6.2).
% SI messages (f.si 1) take SIB1's pdsch-TimeDomainAllocationList when given (siOpt.tdra).
% Table 5.1.2.1.1-2, default A, normal CP: rows {mapping, S, L} for dmrs-TypeA-Position 2 and 3; K0 is 0.
tabA = {'A' 2 12 3 11; 'A' 2 10 3 9; 'A' 2 9 3 8; 'A' 2 7 3 6; 'A' 2 5 3 4; 'B' 9 4 10 4; 'B' 4 4 6 4; ...
        'B' 5 7 5 7; 'B' 5 2 5 2; 'B' 9 2 9 2; 'B' 12 2 12 2; 'A' 1 13 1 13; 'A' 1 6 1 6; 'A' 2 4 2 4; ...
        'B' 4 7 4 7; 'B' 8 4 8 4};
r = tabA(f.tdra + 1, :);
if dmrsPos == 2
    S = r{2}; L = r{3};
else
    S = r{4}; L = r{5};
end
if nargin >= 5 && f.si == 1 && ~isempty(siOpt.tdra)
    row = siOpt.tdra(f.tdra + 1, :);
    r{1} = char('B' - (row(1) == 1));                 % 1: type A, 0: type B
    S = row(2); L = row(3);
end
mcsTab = nrPDSCHMCSTables;
m = mcsTab.QAM64Table(mcsTab.QAM64Table.MCSIndex == f.mcs, :);
pdsch = nrPDSCHConfig('NSizeBWP', nRB, 'NStartBWP', 0, 'Modulation', char(m.Modulation), 'NumLayers', 1, ...
                      'MappingType', r{1}, 'SymbolAllocation', [S L], 'PRBSet', f.rbStart:f.rbStart + f.nRB - 1, ...
                      'VRBToPRBInterleaving', logical(f.vrbInterleaved), 'VRBBundleSize', 2, 'NID', ncellid, ...
                      'RNTI', 65535);
pdsch.DMRS.DMRSConfigurationType = 1;
pdsch.DMRS.DMRSLength = 1;
pdsch.DMRS.DMRSAdditionalPosition = 2;
pdsch.DMRS.DMRSTypeAPosition = dmrsPos;
pdsch.DMRS.NumCDMGroupsWithoutData = 2 - (L == 2 && r{1} == 'B');
pdsch.DMRS.DMRSPortSet = 0;
R = m.TargetCodeRate;
end


function [cw, snrDb, presence, tbs] = sib1_pdsch_llr(carrier, grid, pdsch, R)
% Soft bits of one slot's PDSCH, its SNR, and how well its DM-RS matches (per RB and symbol:
% |sum conj(d) r|^2 / (6 sum |r|^2), about 1/6 for noise, (6 snr + 1) / (6 (snr + 1)) with a PDSCH).
[ind, info] = nrPDSCHIndices(carrier, pdsch);
dInd = nrPDSCHDMRSIndices(carrier, pdsch);
dSym = nrPDSCHDMRS(carrier, pdsch);
[h, nVar] = nrChannelEstimate(carrier, grid, dInd, dSym, 'CDMLengths', pdsch.DMRS.CDMLengths);
[sym, hs] = nrExtractResources(ind, grid, h);
eq = nrEqualizeMMSE(sym, hs, nVar);
snrDb = 10 * log10(mean(abs(hs) .^ 2) / nVar);
r = grid(dInd);
z = reshape(conj(dSym) .* r, 6, []);
e = reshape(abs(r) .^ 2, 6, []);
presence = mean(abs(sum(z, 1)) .^ 2 ./ (6 * sum(e, 1)));
c = nrPDSCHDecode(carrier, pdsch, eq, nVar);
cw = c{1};
tbs = nrTBS(pdsch.Modulation, 1, numel(pdsch.PRBSet), info.NREPerPRB, R);
end


function comb = sib1_combine(tries, carrier, mib, opt, nRB, okFail, siOpt)
% Soft combining over occasions (option 'Combine'; see the help above).
comb = struct('used', 0, 'decoded', {{}}, 'presence', [], 'presenceSfn', [], 'rvRule', 'ceil(3k/2) mod 4, k = (SFN/2) mod 4');
dci = tries([tries.dci] & ~[tries.crcOk] & ~arrayfun(@(t) isempty(t.grid0), tries));
% (a) each slot with its own DCI, grouped by SI indicator and transport block size
if ~isempty(dci)
    key = arrayfun(@(t) t.fields.si * 1e6 + t.tbs, dci);
    for k = unique(key)
        grp = dci(key == k);
        dec = [];
        nCopies = 0;
        for t = grp
            c = carrier; c.NSlot = t.slot; c.NFrame = t.sfn;
            [pdsch, R] = sib1_pdsch_config(t.fields, mib.dmrsTypeAPosition, opt.PCI, nRB, siOpt);
            if t.fields.si == 1 && ~isnan(siOpt.c0crb)
                c.NStartGrid = siOpt.c0crb; pdsch.NStartBWP = siOpt.c0crb;
            end
            [cw, snrDb, ~, tbs] = sib1_pdsch_llr(c, t.grid0, pdsch, R);
            if isempty(dec)
                dec = nrDLSCHDecoder('TargetCodeRate', R, 'TransportBlockLength', tbs, ...
                                     'LDPCDecodingAlgorithm', 'Normalized min-sum', 'MaximumLDPCIterationCount', 30);
            end
            [blk, err] = dec(cw, pdsch.Modulation, 1, t.fields.rv);
            nCopies = nCopies + 1;
            comb.used = comb.used + 1;
            fprintf('combine: + SFN %4d slot %d (own DCI: RB %d+%d, RV %d, SI %d, SNR %+.1f dB): %d copies -> CRC %s\n', ...
                    t.sfn, t.slot, t.fields.rbStart, t.fields.nRB, t.fields.rv, t.fields.si, snrDb, nCopies, okFail{~err + 1});
            if ~err
                hex = lower(reshape(dec2hex(bin2dec(reshape(char(blk(:).' + '0'), 8, []).'), 2).', 1, []));
                comb.decoded{end + 1} = struct('sfn', t.sfn, 'slot', t.slot, 'fields', t.fields, 'tbs', tbs, 'hex', hex, ...
                                               'copies', nCopies);
                reset(dec);
                nCopies = 0;
            end
        end
    end
end
% (b) one allocation for all: also the slots without a decoded DCI whose PDSCH DM-RS is there
allDci = tries([tries.dci]);
if isempty(allDci)
    fprintf('combine: no DCI decoded, so no allocation to assume\n');
    return
end
alloc = arrayfun(@(t) t.fields.rbStart * 1e6 + t.fields.nRB * 1e3 + t.fields.tdra * 10 + t.fields.si, allDci);
if numel(unique(alloc)) > 1
    fprintf('combine: the decoded DCIs use %d different allocations: no allocation assumed for slots without a DCI\n', ...
            numel(unique(alloc)));
    fprintf('combine: %d copies combined, %d SIB/SI decoded\n', comb.used, numel(comb.decoded));
    return
end
ref = find([tries.dci], 1);
f0 = tries(ref).fields;
slot = opt.CombineSlot;
if isempty(slot)
    slot = tries(ref).slot;
end
[pdsch, R] = sib1_pdsch_config(f0, mib.dmrsTypeAPosition, opt.PCI, nRB);
tbs = [];
dec = [];
n = 0;
fprintf('combine: allocation of the DCI in SFN %d slot %d (RB %d+%d, TDRA %d, MCS %d), slot %d of every occasion\n', ...
        tries(ref).sfn, tries(ref).slot, f0.rbStart, f0.nRB, f0.tdra, f0.mcs, slot);
for t = tries
    if t.slot ~= slot || isempty(t.grid0) || t.dci
        continue
    end
    carrier.NSlot = t.slot;
    carrier.NFrame = t.sfn;
    [cw, snrDb, pres, tbsK] = sib1_pdsch_llr(carrier, t.grid0, pdsch, R);
    comb.presence(end + 1) = pres;
    comb.presenceSfn(end + 1) = t.sfn;
    if pres < opt.PresenceThreshold
        continue
    end
    if t.dci
        rv = t.fields.rv;
    else
        rv = mod(ceil(1.5 * mod(floor(t.sfn / 2), 4)), 4);
    end
    if isempty(dec)
        tbs = tbsK;
        dec = nrDLSCHDecoder('TargetCodeRate', R, 'TransportBlockLength', tbs, ...
                             'LDPCDecodingAlgorithm', 'Normalized min-sum', 'MaximumLDPCIterationCount', 30);
    end
    [blk, err] = dec(cw, pdsch.Modulation, 1, rv);
    n = n + 1;
    comb.used = comb.used + 1;
    fprintf('combine: + SFN %4d (RV %d, DM-RS %.2f, SNR %+.1f dB): %d copies -> CRC %s\n', t.sfn, rv, pres, snrDb, n, ...
            okFail{~err + 1});
    if ~err
        hex = lower(reshape(dec2hex(bin2dec(reshape(char(blk(:).' + '0'), 8, []).'), 2).', 1, []));
        comb.decoded{end + 1} = struct('sfn', t.sfn, 'slot', t.slot, 'fields', f0, 'tbs', tbs, 'hex', hex, 'copies', n);
        reset(dec);
        n = 0;
    end
end
fprintf('combine: %d occasions with a PDSCH (of %d in slot %d), %d SIB1 decoded\n', comb.used, numel(comb.presence), ...
        slot, numel(comb.decoded));
end
