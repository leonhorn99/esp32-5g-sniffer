function [rx, meta] = c5_rx_model(tx, fsTx, nRE, opt)
%C5_RX_MODEL  Simulated ESP32-C5 burst capture of a waveform (simulation, not a measurement).
%   [rx, meta] = c5_rx_model(tx, fsTx, nRE, 'SNRdB', -3, 'SampleRate', 20e6, ...)
%   tx    transmitted waveform with SS_REF at 0 Hz (from ssb_test_waveform)
%   fsTx  its sample rate
%   nRE   number of SSB resource elements in tx (all of equal power)
%   rx    N-by-1 complex double in ADC counts (what esp_load_session returns)
%
% Model, in this order:
%   1. resample to the C5 rate (ideal anti-alias filter; the real 4/8 MS/s modes
%      have none, so wideband noise and neighbours fold in: not modelled here)
%   2. frequency shift by FrequencyOffset (SSB centre minus LO) minus the LO error
%      CarrierFrequency*PPM*1e-6 (the C5 crystal error; sample-clock error ignored)
%   3. cut the burst: N samples from StartSample (in C5 samples)
%   4. AWGN at SNRdB, defined per resource element like an SS-SINR:
%      SNR = E_RE / N0 with E_RE the useful-part energy of one SSB RE
%   5. IQ imbalance (gain, phase), scaling to RMSCounts, DC offset,
%      rounding and clipping to Bits (8: CAP16, 10: CAP20)
arguments
    tx (:,1) double
    fsTx (1,1) double
    nRE (1,1) double
    opt.SampleRate (1,1) double = 20e6
    opt.N (1,1) double = 32704
    opt.StartSample double = []          % 1-based, in C5 samples; [] = random
    opt.FrequencyOffset (1,1) double = 3e6
    opt.CarrierFrequency (1,1) double = 2140e6
    opt.PPM (1,1) double = 0
    opt.SNRdB (1,1) double = 10
    opt.Bits (1,1) double = 10
    opt.RMSCounts double = []            % default: 1/8 of full scale (-18 dBFS)
    opt.DCCounts (1,1) double = 0
    opt.IQGainDB (1,1) double = 0
    opt.IQPhaseDeg (1,1) double = 0
    opt.Seed double = []
    opt.EnergyPerRE double = []          % useful energy per RE (sum|tx|^2/fsTx units); set it
                                         % when tx is only a slice of the waveform (nRE unused then)
end
if ~isempty(opt.Seed), rng(opt.Seed); end
fs = opt.SampleRate;
[p, q] = rat(fs/fsTx, 1e-12);
y = resample(tx, p, q);
% SNR definition: energy over the whole waveform, CP part removed (all SSB
% symbols have the normal CP, ratio 2048/2192 at every numerology)
if isempty(opt.EnergyPerRE)
    Ere = sum(abs(y).^2)/fs*(2048/2192)/nRE;
else
    Ere = opt.EnergyPerRE;
end
N0 = Ere/10^(opt.SNRdB/10);
sigma2 = N0*fs;

fShift = opt.FrequencyOffset - opt.CarrierFrequency*opt.PPM*1e-6;
n = (0:numel(y)-1).';
y = y.*exp(1i*2*pi*fShift/fs*n);
N = opt.N;
if isempty(opt.StartSample)
    s0 = randi(numel(y) - N + 1);
else
    s0 = opt.StartSample;
end
y = y(s0:s0+N-1);
y = y + sqrt(sigma2/2)*(randn(N,1) + 1i*randn(N,1));

g = 10^(opt.IQGainDB/20);
ph = opt.IQPhaseDeg*pi/180;
y = real(y) + 1i*g*(imag(y)*cos(ph) - real(y)*sin(ph));

full = 2^(opt.Bits-1);
rmsTarget = opt.RMSCounts;
if isempty(rmsTarget), rmsTarget = full/8; end
scale = rmsTarget/sqrt(mean(abs(y).^2));
y = y*scale + opt.DCCounts*(1+1i)/sqrt(2);
I = min(max(round(real(y)), -full), full-1);
Q = min(max(round(imag(y)), -full+1), full);     % the device sends -Q in [-full, full-1]
rx = I + 1i*Q;

meta.startSample = s0;
meta.tStart = (s0-1)/fs;              % burst start, seconds from the tx waveform start
meta.fShift = fShift;
meta.cfoTrue = -opt.CarrierFrequency*opt.PPM*1e-6;
meta.scale = scale;
meta.sigma2Counts = sigma2*scale^2;
meta.clipFrac = mean(abs([real(y); imag(y)]) >= full - 1);
end
