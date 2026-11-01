function [phi] = evalEigenfunctionKLE(omega_i, a, i, t)
%{
    Evaluates an eigenfunction for an exponential covariance function
    gaussian process. Can handle vector t values

    omega_i = modal frequency
    a = half domain of function (e.g. -a < t < a)
    i = mode number
    t = spatial/temporal variable to evaluate eigenfunction at (t = nx1),
    if t not provided then will return a function of t
   
%}
    fReturnValue = nargin == 4; % e.g. time provided

    commonFactor = sin(2*omega_i.*a) ./ (2*omega_i);
    
    % Faster to duplicate code?
    if fReturnValue
        if isEven(i)
            phi = phiEven(omega_i,t,a,commonFactor);       
        else
            phi = phiOdd(omega_i,t,a,commonFactor);
        end
    else
        % Build functions
        if isEven(i)
            phi = @(t) phiEven(omega_i,t,a,commonFactor);   
        else
            phi = @(t) phiOdd(omega_i,t,a,commonFactor);
        end
    end
end

function f = phiEven(omega_i,t,a,commonFactor)
    f = sin(omega_i.*t) ./ (a - commonFactor).^0.5;
end

function f = phiOdd(omega_i,t,a,commonFactor)
    f = cos(omega_i.*t) ./ (a + commonFactor).^0.5;
end

