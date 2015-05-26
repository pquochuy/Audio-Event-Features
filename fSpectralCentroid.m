function [c] = fSpectralCentroid(x, fs)

% framelength
framelen = length(x);
fftSize = 2^nextpow2(framelen);

X = fft(x,fftSize);
X = abs(X(1:fftSize/2));
X = X';

% frequency tick
f = fs/fftSize * [1 : length(X)];
c = sum(f .* X) / sum(X);
