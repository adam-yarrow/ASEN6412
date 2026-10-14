function [prob] = calculateProbabilityOfMax(muUmax, varUmax, uStats, sigmaLevel)
%{
    Calculates the probability of achieve a maximum value as a parameter of
    the number of MonteCarlo runs.

    Assumes that the cumulative mean and variance for all N runs (e.g. the
    sample statistics) are the population statistics).
%}
sigmaUmax = sqrt(varUmax);
Nmc = numel(varUmax);

% Indicator function
Ic = uStats.max > muUmax(end) + sigmaLevel*sigmaUmax(end);
prob = cumsum(Ic) ./ (1:Nmc);
end