function [z] = melfilterbankenergies(x,fs,M)

% Number of band filters
if isempty(M)
    M = 20;
end

% framelength
framelen = length(x);
% FFT size
fftSize = 2^nextpow2(framelen);

% right FFT, returned size is fftSize/2 + 1
X = rfft(x,fftSize);

% ===== mel-scale filterbank ===
K = fftSize/2 + 1;  % size of each filter
%M = 20; % number of filter
hz2mel = @(hz)(1127*log(1+hz/700)); % Hertz to mel warping function
mel2hz = @(mel)(700*exp(mel/1127)-700); % mel to Hertz warping function
[H, freq ] = trifbank( M, K, [0 fs/2], fs, hz2mel, mel2hz );

% log energy
X = abs(X).^2;
z=H*X;         % multiply x by the power spectrum
z = z';