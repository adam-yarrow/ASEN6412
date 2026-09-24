function [f] = isEven(i)
% Checks if a number is even
f = rem(i, 2) == 0;
end