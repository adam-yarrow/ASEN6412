function L = LagrangePoly1D(y, yInterp, i)
%{
    Returns the interpolated values of fy at yInterp for a univariate y
    y = Nx1
    yInterp = Mx1  
    i = ith poly

    fyInterp = nStates x nTimes x M
%}
M = numel(yInterp);
N = numel(y);

L = ones(M,1);

for k=1:N
    if k ~= i
        L = L .* (yInterp - y(k)) ./ (y(i) - y(k));
    end
end

end