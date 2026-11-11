function [stochasticValues, rawParams, Wcorrected] = MapSCCtoDistributions(...
                    stochasticParam, Y, W)
    %{
        Maps Sparse Smolyak Clenshaw-Curtis weights and points to
        appropriate distributions.

        stochasticParams = structure defining the stochastic parameters
        Y = SCC points on domain Y = [-1,1], nY x nParams
        W = SCC weights

        Outputs:
        stochasticValues = cell array of values for each SCC point
            corresponds to a scaled draw of Y)
        rawParams = cell array of raw values that make up stochasticValues
        Wcorrected = weights corrected by the density of the random
            variable in each dimension at the particular value
    %}

    nParams = numel(stochasticParam);
    nPoints = numel(W);
    rawParams = cell(nParams,1); %% TODO - might not need to be a cell array if working exclusively with numbers?
    
    pdfWeightings = NaN(nPoints, nParams);

    for iParam = 1:nParams
        param = stochasticParam{iParam};
        switch param.distribution
            case 'uniform'
                % Translate Points into domain
                paramMean = (param.min + param.max) / 2;
                paramGradient = (param.max - param.min) / 2;
                rawParams{iParam} = paramMean + paramGradient .* Y(:, iParam);

                % Get Weights
                warning('CHECK IF THIS CONSTANT 1/2 is correct?');
                pdfWeightings(:, iParam) = ones(nPoints,1) * 1/2; % Constant 1/2 factor???

            case 'gaussian'
                error('Not sure if this will work with just gaussian');

            case 'trunc-gaussian'
                % Reference
                % https://en.wikipedia.org/wiki/Truncated_normal_distribution

                % Translate Points into domain
                rawParams{iParam} = param.mean + ...
                    param.truncatedSigmaLevel*param.sigma * Y(:, iParam);
              
                % Get Weights
                pdfWeightings(:, iParam) = param.truncatedSigmaLevel * param.sigma * ...
                    truncNormalPdf(rawParams{iParam}, param.mean, ...
                        param.sigma, param.truncatedSigmaLevel);

         
        end
    end
    
    % Now rescale weights by taking product of all of them along the input
    % dimension direction
    Wcorrected = W .* prod(pdfWeightings,2);   

    % Checking weights are reasonable values given distributions at hand
    % #Use protection it'll save your Stochastic Collocation from spawning
    % terrible children.
    if abs(sum(Wcorrected) - 1) > 0.2
        fprintf(' WARNING: Weights summed to: %f\n',sum(Wcorrected));  
    end  
   
    % Repackage raw data into cell arrays
    stochasticValues = arrayfun(@(i) cellfun(@(x) x(i), rawParams, 'UniformOutput', false), ...
               (1:nPoints)', 'UniformOutput', false);
end