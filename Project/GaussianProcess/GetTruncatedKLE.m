function [d, dataKLE] = GetTruncatedKLE(param)
    %{
        Generates a truncated KLE for a variable of interest (most likely
        Cm_alpha)
    
        Assumes you want to approximate a gaussian process G(t,w), with
        exponential covariance function Gxx = sigma^2*exp(-|t1 - t2|/l)
    %}

    a = (param.inputVariableMax - param.inputVariableMin) / 2;

    [d, lambda, omega] = findMinNumLambdaKLE(param.truncErrorMSE, ...
                                    a, param.l, param.sigma, ...
                                    param.lambdaOptions);

    dataKLE = generateKLE_Struct(lambda, omega, param.meanFunction, a, d);
end