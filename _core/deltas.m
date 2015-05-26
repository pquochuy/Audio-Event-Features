function [d] = deltas(x)
% calculate 1st temporal derivative of matrix x
d = x - [zeros(1,size(x,2)); x(1:end-1,:)];
end