function test_sib1_offline
%TEST_SIB1_OFFLINE  sib1_offline on the 23 Nov band survey: PCI 317 (n7, SSB 2667.75 MHz) at LO 2670.
%   Run from the project root: matlab -batch "addpath('matlab','tests'); test_sib1_offline"
%   Needs data/20261123_165332_survey_2175_2715 (not in samples/: skipped without it).
%   Expected (24 Nov 2026): 5 anchors, CORESET#0 24 RB x 2 symbols at 2665.50-2669.82 MHz, Type0-PDCCH in
%   slots 3/4 of even frames; SIB1 found in slot 3 of SFN 240 and 336 (DCI AL4, RB 0-19, TDRA 1, MCS 4,
%   RV 0, TBS 888), both CRC ok and identical; MCC 232 at bit 42 of the message.
root = fileparts(fileparts(mfilename('fullpath')));
folder = fullfile(root, 'data', '20261123_165332_survey_2175_2715');
if ~isfolder(folder)
    fprintf('test_sib1_offline: SKIPPED (%s not present)\n', folder);
    return
end
out = tempname;
cleanup = onCleanup(@() delete([out '.mat'], [out '.json']));
res = sib1_offline(folder, 'FreqMHz', 2670, 'SSBMHz', 2667.75, 'PCI', 317, 'CFOCenter', -7.7e3, 'Out', out);
assert(numel(res.anchors) == 5 && res.clockResidRms_us < 0.5, 'frame clock from the anchors');
c = res.coreset0;
assert(c.nRB == 24 && c.nSym == 2 && abs(c.fLow - 2665.50e6) < 1, 'CORESET#0 from pdcch-ConfigSIB1 16');
assert(res.type0.n0 == 3 && res.type0.frameParity == 0, 'Type0-PDCCH slots');
ok = res.tries([res.tries.crcOk]);
assert(isequal(sort([ok.sfn]), [240 336]) && all([ok.slot] == 3), 'SIB1 in slot 3 of SFN 240 and 336');
assert(all([ok.tbs] == 888) && strcmp(ok(1).hex, ok(2).hex), 'TBS 888, both copies identical');
f = ok(1).fields;
assert(f.rbStart == 0 && f.nRB == 20 && f.tdra == 1 && f.mcs == 4 && f.rv == 0 && f.si == 0, 'DCI fields');
bits = reshape(dec2bin(hex2dec(cellstr(reshape(ok(1).hex, 2, []).')), 8).', 1, []);
assert(strcmp(bits(43:54), '001000110010'), 'MCC 2 3 2 at bit 42');
fprintf('test_sib1_offline: OK (SIB1 of PCI 317, %d bytes, CRC ok in SFN %s)\n', numel(ok(1).hex) / 2, mat2str([ok.sfn]));
end
