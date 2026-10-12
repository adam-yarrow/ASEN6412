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
    if isEven(i)
        phi = @(t) sin(omega_i.*t) ./ (a - commonFactor).^0.5;       
    else
        phi = @(t) cos(omega_i.*t) ./ (a + commonFactor).^0.5;
    end

    if fReturnValue
        phi = phi(t);
    end
end

