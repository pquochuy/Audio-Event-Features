% count zero crossings
function [z] = fzcr(x)
    x_ = zeros(size(x));
    x_(2 : end) = x(1 : end-1);
    z = sum(abs(sign(x) - sign(x_)))/(2*length(x));
end