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
    omegaGuess = (omegaSolnIdx)*(pi); % skip the first solution as it gives omega = 0
else
    omegaSolnIdx = ceil(i/2);
    fObj = @(w) (1/l) - w*tan(w*a);
    omegaGuess = (omegaSolnIdx - 1)*(pi); 
    % slight perturbation away from zero for i == 1 case 
    if i == 1
        omegaGuess = omegaGuess + 1E-5;
    end
    
end

% Domain bounds for tan function
upperBound = omegaGuess + pi/2;
lowerBound = omegaGuess - pi/2;

% Solve objective functions to find their roots using nonlin lsq as it
% allows for bounding of solution (unlike fsolve)
options = optimoptions('lsqnonlin','Display','none');
[omega, ~, ~, exitFlag] = lsqnonlin(fObj,omegaGuess/a, ...
                            lowerBound/a, upperBound/a, options);
if exitFlag <= 0
    error('Could not solve for omega');
end
end