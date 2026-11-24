function [f, gscn] = ssb_gscn(fmin, fmax)
%SSB_GSCN  NR synchronization raster (GSCN) positions between fmin and fmax (Hz).
%   [f, gscn] = ssb_gscn(fmin, fmax) returns the SS_REF frequencies (Hz, column)
%   and their GSCN numbers, TS 38.101-1 Table 5.4.3.1-1:
%     0-3000 MHz:     SS_REF = N*1200 kHz + M*50 kHz, N = 1..2499, M = 1,3,5,
%                     GSCN = 3N + (M-3)/2
%     3000-24250 MHz: SS_REF = 3000 MHz + N*1.44 MHz, N = 0..14756, GSCN = 7499 + N
%   Band-specific GSCN ranges (e.g. n1: 5279-5419) are not applied here.
%
%   Example: n1 positions the ESP32-C5 can tune to (>= 2130 MHz)
%     [f, g] = ssb_gscn(2130e6, 2170e6);
f = []; gscn = [];
N = (1:2499).';
for M = [1 3 5]
    fm = N*1200e3 + M*50e3;
    f = [f; fm]; %#ok<AGROW>
    gscn = [gscn; 3*N + (M-3)/2]; %#ok<AGROW>
end
N = (0:14756).';
f = [f; 3000e6 + N*1.44e6];
gscn = [gscn; 7499 + N];
keep = f >= fmin & f <= fmax;
[f, order] = sort(f(keep));
g = gscn(keep);
gscn = g(order);
end
