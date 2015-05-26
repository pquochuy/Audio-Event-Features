function [ffbe] = fFFBE(x,fs,M)
% Number of frequency bands
if isempty(M)
    M = 13;
end

z = melfilterbanklogenergies(x,fs,M);

% padding zeros both sides
z = [0,z,0];
% filter with H = z - z^(-1)
h = [-1, 0, 1];

ffbe = zeros(1,M);
for i = 1:M
    ffbe(i) = z(i:i+2)*h';
end
