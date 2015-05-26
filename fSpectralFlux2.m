function [flux] = fSpectralFlux2(energy, preenergy)

flux = 0;
if (isempty(energy) || isempty(preenergy))
    return;
end
% nomalization
len = min([length(energy),length(preenergy)]);
energy = energy / max(energy);
% FFT of previous frame
preenergy = preenergy / max(preenergy);

flux = sum(energy(1:len) - preenergy(1:len));

end