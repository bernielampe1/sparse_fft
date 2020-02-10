function [signalWithNoise, scaledSignal, measuredSNR] = addNoise(signal, SNR)

n = length(signal);
noise = randn(1, n) .* exp(2*pi*1i*rand(1, n));

stdnoise = sqrt(var(noise));
stdsignal = sqrt(var(signal));

scaledSignal = stdnoise / stdsignal * sqrt(10^(SNR/10)) * signal;

signalPower = var(scaledSignal);
noisePower = var(noise);

SNRratio = signalPower / noisePower;
measuredSNR = 10 * log10(SNRratio);
signalWithNoise = scaledSignal + noise';

end
