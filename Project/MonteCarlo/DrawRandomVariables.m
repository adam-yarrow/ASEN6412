function [stochasticValues, rawParams] = DrawRandomVariables(stochasticParams, nSamples)
    %{
        Generates Nsamples of random variables

        stochasticParams = cell array containing structures defining each
        random variables distribution

        stochasticValues = cell array of size Nsamples x 1. Each cell array
        contains a separate array of size nParams x 1 which is each value
        of the parameter of interest for the current sample case

        rawParams = cell array of size Nparam 
    %}

    nParams = numel(stochasticParams);
    rawParams = cell(nParams,1);

    for iParam = 1:nParams
        param = stochasticParams{iParam};
        switch param.distribution
            case 'uniform'
                rawParams{iParam} = param.min + (param.max - param.min)*rand(nSamples,1);
            case 'gaussian'
                rawParams{iParam} = param.mean + param.sigma * randn(nSamples, 1);
            case 'trunc-gaussian'
                x = trandn(-param.truncatedSigmaLevel*ones(nSamples,1), ...
                            param.truncatedSigmaLevel*ones(nSamples,1));
                rawParams{iParam} = param.mean + param.sigma*x;
        end
    end
 
    % Repackage raw data into cell arrays
    stochasticValues = arrayfun(@(i) cellfun(@(x) x(i), rawParams, 'UniformOutput', false), ...
               (1:nSamples)', 'UniformOutput', false);
      
end