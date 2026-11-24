function iq = esp_unpack_iq(payload, n, fmt)
%ESP_UNPACK_IQ  Decode an ESP-SDR CAP16/CAP20 payload to complex samples (ADC counts).
%   iq = esp_unpack_iq(payload, n, fmt)
%   payload  uint8 vector as received (bursts.bin slice)
%   n        number of samples
%   fmt      'CAP16' (I, -Q as int8 pairs) or 'CAP20' (20-bit samples packed LSB
%            first: I = bits 0-9, -Q = bits 10-19, signed 10 bit)
%   iq       n-by-1 complex double, Q already negated (standard orientation)
%
%   Same decoding as capture/esp_sdr.py unpack_iq; tests/test_matlab_loader.m
%   checks both against the fake board's pattern.

payload = uint8(payload(:));
switch fmt
    case 'CAP16'
        need = 2*n;
        assert(numel(payload) >= need, 'esp_unpack_iq:short', 'payload too short');
        v = double(typecast(payload(1:need), 'int8'));
        v = reshape(v, 2, n);
        iq = v(1,:).' - 1i*v(2,:).';
    case 'CAP20'
        need = ceil(n*20/8);
        assert(numel(payload) >= need, 'esp_unpack_iq:short', 'payload too short');
        groups = ceil(n/2);
        b = zeros(5*groups, 1, 'uint64');
        b(1:need) = uint64(payload(1:need));
        b = reshape(b, 5, groups);
        w = b(1,:) + bitshift(b(2,:), 8) + bitshift(b(3,:), 16) + ...
            bitshift(b(4,:), 24) + bitshift(b(5,:), 32);
        s = zeros(2, groups, 'uint64');
        s(1,:) = bitand(w, uint64(hex2dec('FFFFF')));
        s(2,:) = bitand(bitshift(w, -20), uint64(hex2dec('FFFFF')));
        s = s(:);
        s = s(1:n);
        I = sign10(bitand(s, uint64(1023)));
        mQ = sign10(bitand(bitshift(s, -10), uint64(1023)));
        iq = I - 1i*mQ;
    otherwise
        error('esp_unpack_iq:fmt', 'unknown format %s', fmt);
end
end

function v = sign10(u)
v = double(u);
v(v >= 512) = v(v >= 512) - 1024;
end
