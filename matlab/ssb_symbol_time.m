function t = ssb_symbol_time(l, scs)
%SSB_SYMBOL_TIME  Start of the useful part (after the CP) of OFDM symbol l, in seconds.
%   t = ssb_symbol_time(l, scs)
%   l    symbol index counted from the start of the half frame (0-based), any shape
%   scs  subcarrier spacing in kHz (15, 30, ...); normal CP
%   The first symbol of every 0.5 ms has the CP extended by 16*kappa*Tc (TS 38.211 5.3.1).
Tc = 1/(480e3*4096);
kappa = 64;
mu = log2(scs/15);
Nu = 2048*kappa*2^(-mu);
Ncp = 144*kappa*2^(-mu);
Nlong = Ncp + 16*kappa;
per05 = 7*2^mu;                      % symbols per 0.5 ms
h = floor(l/per05);
j = mod(l, per05);
tc = h*round(0.5e-3/Tc) + (j > 0).*((Nlong + Nu) + max(j-1, 0)*(Ncp + Nu));
cp = Ncp + (j == 0)*(Nlong - Ncp);
t = (tc + cp)*Tc;
end
