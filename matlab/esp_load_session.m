function [bursts, meta] = esp_load_session(folder, opts)
%ESP_LOAD_SESSION  Read a capture session written by capture/capture_bursts.py.
%   [bursts, meta] = esp_load_session(folder)
%   [bursts, meta] = esp_load_session(folder, 'Index', 1:10, 'OnlyCrcOk', true)
%   [bursts, meta] = esp_load_session(folder, 'FreqMHz', [2610 2625])   % only bursts at these LOs
%   [recs, meta] = esp_load_session(folder, 'MetaOnly', true)           % records without .iq
%
%   bursts  struct array, one per burst, with the JSON record fields plus
%             .iq     n-by-1 complex double in ADC counts (Q negated)
%             .fs     sample rate in Hz
%             .fLO    LO frequency in Hz (freq_mhz * 1e6)
%   meta    session.json (device answers and settings)
arguments
    folder (1,:) char
    opts.Index double = []
    opts.OnlyCrcOk (1,1) logical = true
    opts.FreqMHz double = []
    opts.MetaOnly (1,1) logical = false
end
meta = jsondecode(fileread(fullfile(folder, 'session.json')));
lines = splitlines(strtrim(fileread(fullfile(folder, 'bursts.jsonl'))));
lines = lines(~cellfun(@isempty, lines));
recs = cellfun(@jsondecode, lines, 'UniformOutput', false);
if opts.OnlyCrcOk
    recs = recs(cellfun(@(r) logical(r.crc_ok), recs));
end
if ~isempty(opts.FreqMHz)          % before loading: a survey session can hold GBs of I/Q
    recs = recs(cellfun(@(r) ismember(round(double(r.freq_mhz)), opts.FreqMHz), recs));
end
if ~isempty(opts.Index)
    recs = recs(opts.Index);
end
fid = fopen(fullfile(folder, 'bursts.bin'), 'r');
cleaner = onCleanup(@() fclose(fid));
bursts = struct([]);
for k = 1:numel(recs)
    r = recs{k};
    if ~opts.MetaOnly
        fseek(fid, r.offset, 'bof');
        payload = fread(fid, r.nbytes, '*uint8');
        r.iq = esp_unpack_iq(payload, r.n, r.fmt);
    end
    r.fLO = double(r.freq_mhz) * 1e6;
    r.fs = double(r.fs);
    r = normalize_empty(r);
    if isempty(bursts)
        bursts = r;
    else
        bursts(end+1) = orderfields(r, bursts(1)); %#ok<AGROW>
    end
end
end

function r = normalize_empty(r)
% JSON null becomes [] in jsondecode; keep the fields so struct arrays line up.
for f = {'t0_us', 'seq', 'lag', 'temp_c', 'gain', 'bandwidth_mhz', 'rms', 'clip_frac'}
    if ~isfield(r, f{1})
        r.(f{1}) = [];
    end
end
end
