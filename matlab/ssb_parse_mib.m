function mib = ssb_parse_mib(trblk, lsbofsfn, hrf, msbidxoffset, Lmax)
%SSB_PARSE_MIB  MIB fields from the outputs of nrBCHDecode (TS 38.331 MIB, TS 38.212 7.1.1).
%   mib = ssb_parse_mib(trblk, lsbofsfn, hrf, msbidxoffset, Lmax)
%   trblk 24 bits: [message choice, SFN 6 MSBs, subCarrierSpacingCommon,
%   ssb-SubcarrierOffset (4 bits), dmrs-TypeA-Position, pdcch-ConfigSIB1 (8 bits),
%   cellBarred, intraFreqReselection, spare]
b = double(trblk(:));
s = double(lsbofsfn(:));
bits2int = @(v) sum(v(:).' .* 2.^(numel(v)-1:-1:0));
mib.messageChoice = b(1);                     % 0 = mib
mib.sfn = bits2int([b(2:7); s]);
mib.scsCommonFlag = b(8);                     % 0: 15 or 60 kHz, 1: 30 or 120 kHz
kssbLsb = bits2int(b(9:12));
if Lmax == 64
    mib.kSSB = kssbLsb;
    mib.ssbIndexMsb = bits2int(msbidxoffset);
else
    mib.kSSB = kssbLsb + 16*double(msbidxoffset(1));
    mib.ssbIndexMsb = [];
end
mib.dmrsTypeAPosition = 2 + b(13);
mib.pdcchConfigSIB1 = bits2int(b(14:21));
mib.cellBarred = b(22) == 0;                  % ENUMERATED {barred, notBarred}
mib.intraFreqReselectionAllowed = b(23) == 0; % ENUMERATED {allowed, notAllowed}
mib.hrf = double(hrf);
end
