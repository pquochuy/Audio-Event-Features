function [sbe,dist] = fSubbandEnergies2(logenergy,nSub)

% Number of sub-bands
if isempty(nSub)
    nSub = 4;
end

M = length(logenergy);

subLen = M/nSub;

sbe = zeros(1,nSub);
for i = 1 : nSub
    sbe(i) = sum(logenergy(((i-1)*subLen + 1) : i*subLen));
end

dist = sbe/sum(sbe);