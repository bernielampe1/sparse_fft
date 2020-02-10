function [SNRDB, ERR] = l1_vs_snr(k, lowSNR, highSNR, snrStep)

% size of data
n = 4096;

% create data set
data_freq = zeros(n, 1);
    
% create k magnitudes at different locations
for j = 1:k
    data_freq(randi([1 4096])) = 1; %randi([1 100]);
end

% create time domain data
data_time = ifft(data_freq) .* n;

% SNRS
snrs = lowSNR:snrStep:highSNR;

% allocate
ERR = zeros(length(snrs), 1);
SNRDB = zeros(length(snrs), 1);

% loop over k
for i = 1:length(snrs)
    % create complex noise
    [signalWithNoise, ~, snr] = addNoise(data_time, snrs(i));

    % run sparse fft
    dft = abs(sfft(signalWithNoise, k));

    dft = dft ./ max(dft);

    % compute l1 error
    e = sum(abs(dft - data_freq)) ./ k;
    ERR(i) = e;

    % SNR
    SNRDB(i) = snr;
end

end
