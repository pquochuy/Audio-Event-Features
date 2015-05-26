function [e] = fShortTimeEnergy(x)
    e = sum(abs(x.^2))/length(x);
end