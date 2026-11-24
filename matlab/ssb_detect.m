function [det, info] = ssb_detect(x, fs, opt)
%SSB_DETECT  Find NR SS/PBCH blocks in one I/Q burst; decode PCI, SSB index, CFO and MIB.
%   [det, info] = ssb_detect(x, fs, 'fLO', 2138e6)
%   x    burst samples (complex, any scale, e.g. ADC counts from esp_load_session)
%   fs   sample rate in Hz
%
%   Options
%     fLO            LO frequency of the capture (Hz)
%     fSSB           candidate SSB centres (SS_REF, Hz). Default: every GSCN position
%                    whose whole SSB fits inside +-fs/2 around fLO. fLO and fSSB are
%                    absolute RF frequencies: the symbol phase below depends on them
%     SymbolPhase    undo the phase that TS 38.211 5.4 gives every OFDM symbol (default
%                    true). Each subcarrier's phase restarts at the useful start of every
%                    symbol at its absolute RF frequency, so a receiver tuned to SS_REF
%                    sees a step of 2*pi*fSSB*Tsym from symbol to symbol. Left in, the
%                    PSS->SSS refinement reads it as a CFO of -fSSB mod 7.0 kHz (+3.93 kHz
%                    at 2125.45 MHz, the bias of the first-light results) and the DM-RS
%                    search for i_SSB loses its coherence across symbols
%     SubcarrierSpacing  15 (default, n1/n7 FDD, case A) or 30 (case B/C)
%     BlockPattern   'Case A' (15 kHz default) | 'Case B' | 'Case C'
%     Lmax           SSB candidates per half frame: 4 (f <= 3 GHz) or 8
%     CFOMax, CFOStep  coarse CFO search (Hz). The C5 crystal is specified at
%                    +-10 ppm = +-21 kHz at 2.14 GHz
%     CFOCenter      centre of that search (Hz), e.g. the board's known crystal error
%                    times the frequency (at 4.9 GHz it can exceed 30 kHz)
%     PSSThreshold   PSS statistic threshold. The statistic is unit-mean exponential
%                    under white noise for one (lag, CFO, N_ID2) trial
%     SSSThreshold   SSS statistic threshold (same normalisation, per N_ID1 trial)
%     DecodeMIB, Combine, ListLength   PBCH/BCH decoding and soft combining of the
%                    SSBs of one half frame inside the burst
%     t0_us          burst start (microseconds since boot, from CAPT); adds absolute times
%
%   det  struct array, one row per PSS peak that passed PSSThreshold:
%        fSSB, nid2, pssStat, cfo, n0 (1-based sample in x, fractional), tPss (s from
%        x(1)), full (whole SSB inside the burst), sssStat, sssRatio, sssOk, pci, ibar,
%        iSSB, dmrsMetric, snrDb, mibOk, mib, tHalfFrame, tFrame, (tFrame_us), sssAll (the
%        336 SSS statistics of a full SSB, for accumulation over bursts in detect_session)
%   info candidate list, combined decodes (info.combined)
%
%   Method: per candidate, mix to 0 Hz and resample to 512*SCS; PSS search over CFO
%   hypotheses on a +-1.1 MHz (x SCS/15) low-passed copy; fine CFO from the two PSS
%   halves; SSS (coherent with the PSS channel) -> PCI; CFO refinement from the
%   PSS->SSS phase (after undoing the 38.211 symbol phase); PBCH DM-RS -> i_SSB (and half frame for Lmax 4); MMSE equalisation
%   and polar decoding of the PBCH/BCH (5G Toolbox).
arguments
    x (:,1) double
    fs (1,1) double
    opt.fLO (1,1) double = 0
    opt.fSSB double = []
    opt.SubcarrierSpacing (1,1) double = 15
    opt.BlockPattern (1,:) char = ''
    opt.Lmax (1,1) double = 4
    opt.CFOMax (1,1) double = 30e3
    opt.CFOStep (1,1) double = 5e3
    opt.CFOCenter (1,1) double = 0
    opt.PSSThreshold (1,1) double = 18
    opt.SSSThreshold (1,1) double = 13
    opt.MaxPeaks (1,1) double = 6
    opt.SmoothSC (1,1) double = 9        % subcarriers of channel smoothing for SSS/DM-RS (135 kHz)
    opt.DecodeMIB (1,1) logical = true
    opt.Combine (1,1) logical = true
    opt.ListLength (1,1) double = 8
    opt.t0_us double = []
    opt.KeepTrace (1,1) logical = false  % info.trace: PSS statistic over time per candidate (figures)
    opt.SymbolPhase (1,1) logical = true
end
scs = opt.SubcarrierSpacing;
pattern = opt.BlockPattern;
if isempty(pattern)
    pattern = 'Case A';
    if scs == 30, pattern = 'Case C'; end
end
symFirst = first_symbols(pattern, opt.Lmax);
Nfft = 512;
fsS = Nfft*scs*1e3;
cpLen = 36;                        % normal CP at Nfft 512 (all SSB symbols)
symLen = Nfft + cpLen;
backoff = 9;                       % FFT window starts this many samples into the CP
halfBW = 120*scs*1e3;

x = x - mean(x);                   % static DC offset of the ESP receiver
N = numel(x);
if isempty(opt.fSSB)
    half = fs/2 - halfBW - 0.1e6;
    if half <= 0
        error('ssb_detect:narrow', 'fs = %g is too narrow for a GSCN search; pass fSSB', fs);
    end
    fC = ssb_gscn(opt.fLO - half, opt.fLO + half);
else
    fC = opt.fSSB(:);
end
[p, q] = rat(fsS/fs, 1e-12);
nx = (0:N-1).';

K = constants(scs, fsS, Nfft);
cfoH = opt.CFOCenter + (-opt.CFOMax:opt.CFOStep:opt.CFOMax);

det = empty_det();
det = det([]);                     % 0x0 with all fields, so [det.x] works when nothing is found
traces = struct('fSSB', {}, 'fsS', {}, 'stat', {}, 'cfo', {});
for ic = 1:numel(fC)
    y = x.*exp(-1i*2*pi*(fC(ic) - opt.fLO)/fs*nx);
    z = resample(y, p, q);
    L = numel(z);
    if L < Nfft + 2*symLen, continue; end
    zp = conv(z, K.hLP, 'same');
    E = movsum(abs(zp).^2, [0 Nfft-1]);
    nLag = L - Nfft + 1;
    E = E(1:nLag);
    Nf = 2^nextpow2(L + Nfft);
    Rf = conj(fft(K.r, Nf));
    tz = (0:L-1).';
    best = zeros(nLag, 3);
    bestH = ones(nLag, 3);
    for ih = 1:numel(cfoH)
        C = ifft(fft(zp.*exp(-1i*2*pi*cfoH(ih)/fsS*tz), Nf).*Rf);
        S = abs(C(1:nLag, :)).^2 ./ ((E/Nfft).*K.rE) * K.enbwFrac;
        upd = S > best;
        best(upd) = S(upd);
        bestH(upd) = ih;
    end
    if opt.KeepTrace
        traces(end+1) = struct('fSSB', fC(ic), 'fsS', fsS, 'stat', best, 'cfo', cfoH(bestH)); %#ok<AGROW>
    end
    for nid2 = 0:2
        s = best(:, nid2+1);
        for pk = 1:opt.MaxPeaks
            [m, i0] = max(s);
            if m < opt.PSSThreshold, break; end
            s(max(1, i0-symLen):min(nLag, i0+symLen)) = 0;
            d = process_peak(z, zp, i0, nid2, m, cfoH(bestH(i0, nid2+1)), K, opt, ...
                fsS, Nfft, symLen, backoff, symFirst, scs, opt.SymbolPhase*fC(ic));
            d.fSSB = fC(ic);
            d.n0 = (i0 - 1)*fs/fsS + 1;
            if ~isempty(opt.t0_us)
                d.tPss_us = opt.t0_us + d.tPss*1e6;
                d.tFrame_us = opt.t0_us + d.tFrame*1e6;
            end
            det = append_row(det, d);
        end
    end
end
det = dedupe(det, 3*fs/fsS);

info.candidates = fC;
info.trace = traces;
info.fsSearch = fsS;
info.combined = struct([]);
if opt.Combine && opt.DecodeMIB && ~isempty(det)
    [det, info.combined] = combine_half_frames(det, opt);
end
det = rmfield(det, 'llr');
end

% -------------------------------------------------------------------------
function d = process_peak(z, zp, n0, nid2, stat, cfoC, K, opt, fsS, Nfft, symLen, b, symFirst, scs, fRF)
d = empty_det();
d.nid2 = nid2;
d.pssStat = stat;
d.tPss = (n0 - 1)/fsS;
r = K.r(:, nid2+1);
t = (n0-1:n0+Nfft-2).';
seg = zp(n0:n0+Nfft-1).*exp(-1i*2*pi*cfoC/fsS*t);
h = Nfft/2;
c1 = sum(seg(1:h).*conj(r(1:h)));
c2 = sum(seg(h+1:end).*conj(r(h+1:end)));
cfo = cfoC + angle(c2*conj(c1))/(2*pi*h/fsS);
d.cfo = cfo;
L = numel(z);
d.full = (n0 - b >= 1) && (n0 + 3*symLen + Nfft - 1 - b <= L);
if ~d.full
    return
end
G = demod(z, n0, cfo, fsS, Nfft, symLen, b, fRF);
% SSS, coherent with the channel seen by the PSS (same 127 subcarriers). The PSS
% estimate is smoothed over SmoothSC subcarriers: without it the product of two noisy
% terms costs about 5 dB at -5 dB SNR (simulation of 19 Nov, first version).
H = movmean(G(57:183, 1).*nrPSS(nid2), opt.SmoothSC);
a = G(57:183, 3).*conj(H);
cS = K.S{nid2+1}.'*a;
st = abs(cS).^2/sum(abs(a).^2);
d.sssAll = st;                      % all 336 N_ID1 statistics, for accumulation over bursts
[st, order] = sort(st, 'descend');
d.sssStat = st(1);
d.sssRatio = st(1)/st(2);
d.sssOk = st(1) >= opt.SSSThreshold;
if ~d.sssOk
    return
end
nid1 = order(1) - 1;
d.pci = 3*nid1 + nid2;
% CFO refinement: PSS -> SSS phase over two symbols. It is unambiguous only within
% +-fsS/(4 symLen) = +-3.5 kHz, and the PSS-halves estimate before it can be off by more.
% On the first-light captures (2026-11-23) the PBCH then failed with the CFO exactly one
% period (7.0 kHz) off; most of that came from the 38.211 symbol phase, which demod now
% removes (fixed 23 Nov afternoon). The refined value is used first; only if its MIB fails
% the CRC are the two candidates one period away tried, and one of them is kept only if
% its MIB passes the CRC. (Choosing by the DM-RS metric alone picked wrong candidates in
% the AWGN simulation, 23 Nov.)
cfo0 = cfo + angle(cS(order(1)))/(2*pi*2*symLen/fsS);
period = fsS/(2*symLen);
c = eval_candidate(z, n0, cfo0, d.pci, opt, fsS, Nfft, symLen, b, fRF);
if opt.DecodeMIB && ~c.mibOk
    for kc = [-1 1]
        ck = eval_candidate(z, n0, cfo0 + kc*period, d.pci, opt, fsS, Nfft, symLen, b, fRF);
        if ck.mibOk
            c = ck;
            break
        end
    end
end
d.cfo = c.cfo;
d.ibar = c.ibar;
d.dmrsMetric = c.dmrsMetric;
if opt.Lmax == 4
    d.iSSB = mod(d.ibar, 4);
    hrfDmrs = floor(d.ibar/4);
else
    d.iSSB = d.ibar;
    hrfDmrs = NaN;
end
d.hrfDmrs = hrfDmrs;
d.tHalfFrame = d.tPss - ssb_symbol_time(symFirst(d.iSSB+1), scs);
d.tFrame = d.tHalfFrame - 5e-3*hrfDmrs;
if ~opt.DecodeMIB
    return
end
d.llr = c.llr;
d.snrDb = c.snrDb;
d.mibOk = c.mibOk;
d.mib = c.mib;
d.pbchEq = c.eq;
if d.mibOk
    d.tFrame = d.tHalfFrame - 5e-3*d.mib.hrf;
end
end

function c = eval_candidate(z, n0, cfo, pci, opt, fsS, Nfft, symLen, b, fRF)
% Demodulate the SSB with one CFO value, find i_SSB from the PBCH DM-RS and decode the PBCH.
c = struct('cfo', cfo, 'ibar', NaN, 'dmrsMetric', NaN, 'llr', [], 'snrDb', NaN, 'mibOk', false, 'mib', [], 'eq', []);
G = demod(z, n0, cfo, fsS, Nfft, symLen, b, fRF);
% i_SSB (and half frame for Lmax 4) from the PBCH DM-RS, coherently against the smoothed
% PSS+SSS channel on the subcarriers they share (56..182, symbols 1 and 3)
dmrsInd = double(nrPBCHDMRSIndices(pci));
[scOf, symOf] = ind2sub([240 4], dmrsInd);
use = scOf >= 57 & scOf <= 183 & (symOf == 2 | symOf == 4);
Href = zeros(240, 1);
Href(57:183) = movmean((G(57:183, 1).*nrPSS(pci) + G(57:183, 3).*nrSSS(pci))/2, opt.SmoothSC);
metric = zeros(8, 1);
for ibar = 0:7
    ref = nrPBCHDMRS(pci, ibar);
    v = G(dmrsInd(use)).*conj(ref(use)).*conj(Href(scOf(use)));
    metric(ibar+1) = abs(sum(v))^2/sum(abs(v).^2);
end
[c.dmrsMetric, k] = max(metric);
c.ibar = k - 1;
if ~opt.DecodeMIB
    return
end
refInd = [dmrsInd; nrSSSIndices; nrPSSIndices];
refSym = [nrPBCHDMRS(pci, c.ibar); nrSSS(pci); nrPSS(pci)];
% nrChannelEstimate wants a whole slot without a carrier object: pad to 14 symbols
% (column-major, so the 240x4 linear indices stay valid); only pilots are used.
G14 = zeros(240, 14);
G14(:, 1:4) = G;
[hest, nVar] = nrChannelEstimate(G14, refInd, refSym);
hest = hest(:, 1:4);
pbchInd = nrPBCHIndices(pci);
[eq, csi] = nrEqualizeMMSE(G(pbchInd), hest(pbchInd), nVar);
c.eq = eq;                          % equalised PBCH symbols (for constellation figures)
if opt.Lmax == 4
    v = mod(c.ibar, 4);
else
    v = c.ibar;
end
c.llr = nrPBCHDecode(eq, pci, v, nVar).*repelem(csi, 2);
c.snrDb = 10*log10(mean(abs(hest(pbchInd)).^2)/nVar);
[c.mibOk, c.mib] = decode_bch(c.llr, opt, pci);
end

function [ok, mib] = decode_bch(llr, opt, pci)
[~, err, trblk, lsbsfn, hrf, msb] = nrBCHDecode(llr, opt.ListLength, opt.Lmax, pci);
ok = ~err;
mib = [];
if ok
    mib = ssb_parse_mib(trblk, lsbsfn, hrf, msb, opt.Lmax);
    mib.trblk = trblk;               % decoded 24 MIB bits (for re-encoding checks)
end
end

function G = demod(z, n0, cfo, fsS, Nfft, symLen, b, fRF)
% The four SSB symbols, CFO removed, FFT window b samples into the CP. fRF (absolute
% SS_REF frequency, 0 = off) undoes the TS 38.211 5.4 symbol phase: all four SSB symbols
% have the normal CP, so it is a step of 2*pi*fRF*symLen/fsS per symbol.
k = (0:239).' - 120;
bins = mod(k, Nfft) + 1;
comp = exp(1i*2*pi*k*b/Nfft);      % undo the phase ramp of the early FFT window
G = zeros(240, 4);
for s = 0:3
    idx = (n0 + s*symLen - b : n0 + s*symLen - b + Nfft - 1).';
    F = fft(z(idx).*exp(-1i*2*pi*cfo/fsS*(idx - 1)));
    G(:, s+1) = F(bins).*comp*exp(1i*2*pi*mod(fRF*s*symLen/fsS, 1));
end
end

function [det, comb] = combine_half_frames(det, opt)
% Soft-combine the PBCH LLRs of SSBs of one cell in one half frame (same BCH payload).
comb = struct([]);
use = find(arrayfun(@(d) d.full && d.sssOk && ~isempty(d.llr), det));
done = false(size(det));
for a = use(:).'
    if done(a), continue; end
    grp = use(arrayfun(@(k) det(k).pci == det(a).pci && ...
        abs(det(k).tHalfFrame - det(a).tHalfFrame) < 20e-6 && ...
        isequaln(det(k).hrfDmrs, det(a).hrfDmrs), use));
    done(grp) = true;
    if numel(grp) < 2, continue; end
    llr = sum([det(grp).llr], 2);
    [ok, mib] = decode_bch(llr, opt, det(a).pci);
    c.pci = det(a).pci;
    c.members = grp;
    c.iSSB = [det(grp).iSSB];
    c.mibOk = ok;
    c.mib = mib;
    c.tHalfFrame = mean([det(grp).tHalfFrame]);
    comb = append_row(comb, c);
    for k = grp
        det(k).mibOkComb = ok;
    end
end
end

function det = dedupe(det, tolSamples)
% The same SSB found from neighbouring raster candidates: keep the strongest.
if numel(det) < 2, return; end
[~, order] = sort([det.pssStat], 'descend');
det = det(order);
keep = true(size(det));
for a = 1:numel(det)
    if ~keep(a), continue; end
    for b = a+1:numel(det)
        if keep(b) && abs(det(b).n0 - det(a).n0) <= tolSamples && ...
                (isequaln(det(b).pci, det(a).pci) || isnan(det(b).pci) || det(b).nid2 == det(a).nid2)
            keep(b) = false;
        end
    end
end
det = det(keep);
[~, order] = sort([det.n0]);
det = det(order);
end

function d = empty_det()
d = struct('fSSB', NaN, 'nid2', NaN, 'pssStat', NaN, 'cfo', NaN, 'n0', NaN, 'tPss', NaN, ...
    'full', false, 'sssStat', NaN, 'sssRatio', NaN, 'sssOk', false, 'pci', NaN, ...
    'ibar', NaN, 'iSSB', NaN, 'hrfDmrs', NaN, 'dmrsMetric', NaN, 'snrDb', NaN, ...
    'mibOk', false, 'mibOkComb', false, 'mib', [], 'tHalfFrame', NaN, 'tFrame', NaN, ...
    'tPss_us', NaN, 'tFrame_us', NaN, 'sssAll', [], 'pbchEq', [], 'llr', []);
end

function s = append_row(s, r)
if isempty(s)
    s = r;
else
    s(end+1) = r;
end
end

function K = constants(scs, fsS, Nfft)
persistent cache
key = sprintf('%g_%g', scs, Nfft);
if isempty(cache), cache = containers.Map; end
if isKey(cache, key)
    K = cache(key);
    return
end
k = (56:182).' - 120;               % PSS subcarriers relative to SS_REF
K.r = zeros(Nfft, 3);
for nid2 = 0:2
    X = zeros(Nfft, 1);
    X(mod(k, Nfft) + 1) = nrPSS(nid2);
    K.r(:, nid2+1) = ifft(X);
end
K.rE = sum(abs(K.r).^2, 1);
K.hLP = fir1(128, (1.1e6*scs/15)/(fsS/2)).';
K.enbwFrac = sum(K.hLP.^2);         % filtered-noise power / input-noise power (white input)
K.S = cell(1, 3);
for nid2 = 0:2
    S = zeros(127, 336);
    for nid1 = 0:335
        S(:, nid1+1) = nrSSS(3*nid1 + nid2);
    end
    K.S{nid2+1} = S;
end
cache(key) = K;
end

function s = first_symbols(pattern, L)
switch pattern
    case {'Case A', 'Case C'}
        s = reshape([2; 8] + 14*(0:(L/2 - 1)), 1, []);
    case 'Case B'
        s = reshape([4; 8; 16; 20] + 28*(0:(L/4 - 1)), 1, []);
    otherwise
        error('ssb_detect:pattern', 'pattern %s not supported', pattern);
end
end
