function [lambda] = getLambdaKLE(l, omega)
%{
    Returns the ith eigenvalue of the KL expansion for an exponential
    covariance function of a random process.

    l = length scale
    omega = ith frequency
%}
lambda = 2*l./(l.^2 .* omega.^2 + 1);
end