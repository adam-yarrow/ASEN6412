function [lambda] = getLambdaKLE(l, omega, sigma)
%{
    Returns the ith eigenvalue of the KL expansion for an exponential
    covariance function of a random process.

    l = length scale
    omega = ith frequency
    sigma = standard deviation of the exponential covariance function.
%}
lambda = 2*l./(l.^2 .* omega.^2 + 1)*sigma^2;
end