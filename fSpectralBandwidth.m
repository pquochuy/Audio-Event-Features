function [sb] = fSpectralBandwidth(x, fs)

% framelength
framelen = length(x);
fftSize = 2^nextpow2(framelen);

X = fft(x,fftSize);
X = abs(X(1:fftSize/2));
X = X';

% frequency tick
f = fs/fftSize * [1 : length(X)];

% computing spectral centroid
centroid = sum(f .* X) / sum(X);

sb = sqrt( sum(((f - centroid).^2) .* (X.^2)) / sum(X.^2));