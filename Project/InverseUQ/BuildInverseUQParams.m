function const = BuildInverseUQParams(mode, covInflationSF)

    if nargin < 2
        covInflationSF = 1;
    end
    
    % Baseline params
    const = ModelParams();

    % Modifications to baseline
    switch lower(mode)
        case lower('static-cmAlpha')
            % do nothing, assumes default behaviour
        case lower('linear-cmAlpha')
            const.aero.movMeanCM_alpha.fActive = true;

        case lower('poly-cmAlpha')
            const.aero.polyCM_alpha.fActive = true;

        case lower('Cmalpha-GP-staticMean') 
            const = buildBaselineGPparam(const, covInflationSF);
                
        case lower('Cmalpha-GP-movMean')
            const = buildBaselineGPparam(const, covInflationSF);
            const.aero.movMeanCM_alpha.fActive = true;
    end
    
end

function const = buildBaselineGPparam(const, covInflationSF)
    % Get GP parameters
    [~, stochasticParamGP, ~] = ...
        ForwardUqPerturbationVars('Cmalpha-GP-staticMean', const);
    stochasticParamGP.sigma = stochasticParamGP.sigma * sqrt(covInflationSF);

    % Get Truncated KLE
    [d, dataKLE] = GetTruncatedKLE(stochasticParamGP);

    % Update const struct
    const.aero.GP.fActive = true;
    const.aero.GP.dataKLE = dataKLE;
    Ktrunc = stochasticParamGP.YtruncSigma;
    const.aero.GP.Y = Ktrunc*rand(d,1);
end
