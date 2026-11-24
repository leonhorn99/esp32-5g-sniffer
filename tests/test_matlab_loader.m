function test_matlab_loader
%TEST_MATLAB_LOADER  esp_load_session against the fixture written by tests/fake_device.py.
%   Run from the project root: matlab -batch "addpath('matlab','tests'); test_matlab_loader"
root = fileparts(fileparts(mfilename('fullpath')));
[b, meta] = esp_load_session(fullfile(root, 'tests', 'fixtures', 'session_pattern'));
assert(strcmp(meta.device.info, 'C5SDR 6 burst 16380'));
assert(numel(b) == 2);
% burst 1: CAPT20, n = 101, seq 0; burst 2: CAPT16, n = 64, seq 1 (pattern() in fake_device.py)
specs = {'CAP20', 101, 0, 1024; 'CAP16', 64, 1, 256};
for k = 1:2
    [fmt, n, seq, m] = specs{k, :};
    assert(strcmp(b(k).fmt, fmt) && b(k).n == n && numel(b(k).iq) == n);
    kk = (0:n-1).';
    I = mod(37*kk + 11*seq, m) - m/2;
    Q = mod(91*kk + 5*seq, m - 1) - (m/2 - 1);
    assert(isequal(real(b(k).iq), I), 'I mismatch in burst %d', k);
    assert(isequal(imag(b(k).iq), Q), 'Q mismatch in burst %d', k);
end
assert(b(1).fs == 20e6 && b(2).fs == 4e6 && b(1).fLO == 2138e6);
assert(b(1).t0_us == 1000000 && b(2).t0_us == 1080000);
fprintf('test_matlab_loader: OK\n');
end
