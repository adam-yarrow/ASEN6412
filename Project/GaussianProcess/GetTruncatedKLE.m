function [d, dataKLE] = GetTruncatedKLE(param)
    %{
        Generates a truncated KLE for a variable of interest (most likely
        Cm_alpha)
    
        Assumes you want to approximate a gaussian process G(t,w), with
        exponential covariance function Gxx = sigma^2*exp(-|t1 - t2|/l)
    %}

    a = (param.inputVariableMax - param.inputVariableMin) / 2;

    % [d, lambda, omega] = findMinNumLambdaKLE(param.truncErrorMSE, ...
    %                                 a, param.l, param.sigma, ...
    %                                 param.lambdaOptions);

    % Fixed number of modes
    d = param.d;
    omega = NaN(d,1);
    lambda = NaN(d,1);
    for i = 1:param.d
        omega(i) = getOmegaKLE(i,a,param.l);
        lambda(i) = getLambdaKLE(param.l,omega(i),param.sigma);
    end

    dataKLE = generateKLE_Struct(lambda, omega, param.meanFunction, a, d);
end