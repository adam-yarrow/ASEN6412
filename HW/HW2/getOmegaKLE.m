function [omega] = getOmegaKLE(i, a, l)
%{
    This function returns the frequency for the ith eigenfunction of the
    KL expansion of a process with an exponential covariance matrix.

    i = index for eigenfunction
    a = bounds on domain
    l = length scale
%}

% Check Even/Odd and what solution index we are up to
fIsEven = isEven(i);

if fIsEven  
    omegaSolnIdx  = i/2;
    fObj = @(w) w + (1/l)*tan(w*a);
else
    omegaSolnIdx = ceil(i/2);
    fObj = @(w) (1/l) - w*tan(w*a);
end

% Get initial Guess Starting Location
% FIXME: Maybe change this based on odd or even solution?
omegaGuess = (omegaSolnIdx - 1)*(2*pi) + pi/4; % starting to the right of the middle of the domain.

% Solve objective functions to find their roots using newton-rhapson
options = optimoptions('fsolve','Display','none');
[omega, ~, exitFlag] = fsolve(fObj,omegaGuess,options);

if exitFlag <= 0
    error('Could not solve for omega');
end
end