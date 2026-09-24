function [d, lambda, omega] = findMinNumLambdaKLE(mseThreshold, a, l, sigma, options)
%{
    Calculates the minimum number of eigenfunctions required to obtain an
    approximation of a random process which achieves a specified Mean
    Square Error threshold. 

    Assumes an exponentially decaying covariance function.

    Will return the minimum number of lambda and their values (+omega
    values);

    mseThreshold = threshold as a fraction for acceptable error
    a = domain bound
    l = characteristic length of exp cov function
    sigma = std dev for exp cov function
    options = struct containing
        - mode = 'analytical','approximate'
        - k = (Optional) If approximate, k is number of modes in advance to estimate
        with
        - MSE = (Optional) If analytical this is the analytical MSE.
%}

alpha = 1 - mseThreshold;
truncErrorFrac = 0;
i = 1;
d = 0;

lambda = [];
omega = [];

lambdaSum = 0;

mode = options.mode;

% Iterate until meet threshold
while truncErrorFrac < alpha

    % Analytical Approach
    if strcmpi('analytical',mode)
        omega_i = getOmegaKLE(i,a,l);
        lambda_i = getLambdaKLE(l,omega_i,sigma);
        lambdaSum = lambdaSum + lambda_i;
    
        truncErrorFrac = lambdaSum / options.MSE;
           
        lambda = [lambda, lambda_i];
        omega = [omega, omega_i];
        d = i;
    else

    end
    

    i = i+1;

end


 


end