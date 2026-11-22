function [xStats, alphaStats, rawData] = StochasticColloc(const, Z, W, ...
                        controller, stochasticParam, stochasticParamGP, ...
                        overrideSimData)
    %{
        This function runs stochastic collocation to estimate the
        trajectory statistics and forward UQ.

        There's very little idiot proofing on these inputs, we raw dogging
        a bunch of shit, hence sometimes get some "fun" behaviours.

        Inputs:
            const = baseline parameter set
            Z = d x N points to evaluate at (Z = [-1,1])
            W = N x 1 weights for each point
            controller = annonymous function representing the controller
            stochasticParam = cell array of structs defining how to perturb the
            parameters in the model.
            stochasticParamGP = struct containing gaussian process data for
            Cm_alpha
            overrideSimData = cell array of precomputed sim data (e.g. from
            other levels, will negate the need to run the simulation)

        Ouputs:
            xStats = statistics struct
            alphaStats = like xStats but alpha
            rawData = raw sim data/function evals - handy for nested 
            behavior ("very nice!" - Borat)   
    %}
    % No override sim data case
    if nargin < 7
        fNoExternalSimData = true;
    else
        fNoExternalSimData = false;
    end

    % General
    nAbsisscas = numel(W);

    % Storage
    nTimes = const.tEnd/const.dT + 1;
    xStorage = NaN(const.nStates, nTimes);
    uStorage = NaN(nTimes,1);
    alphaStorage = NaN(nTimes, 1);

    % Paths to perturb
    pathsToPerturb = cellfun(@(x) x.path, stochasticParam, 'UniformOutput', false);

    % Handle Gaussian Processes First
    ZnonGP = Z;
    dGP = 0;
    Ygp = NaN(nAbsisscas, dGP);
    fGP = false;

    if ~isempty(stochasticParamGP)
        % Set general parameters that will be active for all cases
        const.aero.GP.fActive = true;
        const.aero.GP.dataKLE = stochasticParamGP.dataKLE;
        dGP = stochasticParamGP.dataKLE.d;
        fGP = true;

        % Map Zi = [-1,1] to Ygp_i ~ [-K,k]
        Ygp = stochasticParamGP.YtruncSigma * Z(:,1:dGP); 

        % Adjust weights based on truncated KLE - same process as truncated
        % gaussian
        Wgp = NaN(nAbsisscas, dGP);        
        for iDim = 1:dGP
            % Weights assume uniform random variable inputs
            Wgp(:,iDim) = 1/2* ones(size(Ygp(:,iDim)));            
        end

        % Overite old W with new GP corrected W scaled by all dimensions
        W = W .* prod(Wgp, 2);

        % Disabled because full weight product not realized yet
        % % Check W is well conditioned (should be close to 1, checking +/- 20%)
        % if abs(sum(W) - 1) > 0.2
        %     fprintf('WARNING: GP weights summed to: %f\n',sum(W));   
        % end 

        % Remove Z values that were consumed by the dimensions of the
        % gaussian process
        ZnonGP = Z(:,dGP+1:end);
    end

    % Make Stochastic Colloc Inputs
    [stochasticValues, rawStochasticValues, Wcorrected] = MapSCCtoDistributions(...
                    stochasticParam, ZnonGP, W);

    % Checking weights are reasonable values given distributions at hand
    % (Use protection it'll save your Stochastic Collocation from spawning
    % terrible children.)
    if abs(sum(Wcorrected) - 1) > 0.2
        fprintf(' WARNING: Weights summed to: %f\n',sum(Wcorrected));  
    end  
   
    % Run Stoch Colloc
    t0 = tic;
    wb = progressBar(1, nAbsisscas, t0, [], 'SC Runs');
    simDataStorage = cell(nAbsisscas, 1);

    xStats = struct();
    xStats.mean = zeros(const.nStates, nTimes);
    xStats.cov = zeros(const.nStates, const.nStates, nTimes);
    xStats.corr = zeros(const.nStates, const.nStates, nTimes);

    ExSq = zeros(const.nStates, const.nStates, nTimes);

    alphaStats.mean = zeros(nTimes,1);
    alphaStats.var = zeros(nTimes,1);

    EalphaSq = zeros(nTimes, 1);

    % Process each absisscas
    for i = 1:nAbsisscas
        % Overwrite perturbations
        constTemp = updateStruct(const, pathsToPerturb, stochasticValues{i});
            
        % Process GP variables
        if fGP
            constTemp.aero.GP.Y = Ygp(i,:)';            
        end

        % Run Model
        if fNoExternalSimData
            simData = Simulation(controller, constTemp);
        else
            simData = overrideSimData{i};
        end

        simDataStorage{i} = simData;
        xStorage(:,:,i) = simData.x;
        uStorage(:,i) = simData.u;
        alphaStorage(:,i) = calcAlpha(simData.x);

        % Calculate Statistics Using Stochastic Collocation
        % E[x]
        xStats.mean = xStats.mean + simData.x * Wcorrected(i); 

        % E[alpha]
        alphaStats.mean = alphaStats.mean + alphaStorage(:,i)*Wcorrected(i);
        
        % E[x^2]
        % FIXME: Maybe can do with pagemtimes?
        for timeIdx = 1:nTimes
            ExSq(:,:,timeIdx) = ExSq(:,:,timeIdx) + ...
                Wcorrected(i)*simData.x(:,timeIdx)*simData.x(:,timeIdx)';
        end

        % E[alpha^2]
        EalphaSq = EalphaSq + Wcorrected(i)*alphaStorage(:,i).^2;
          
        % Progress bar
        wb = progressBar(i, nAbsisscas, t0, wb);    
    end

    % Process cov and correlation - FIXME: Find a way to get around this
    % for loop
    for timeIdx = 1:nTimes
        % Cov
        xStats.cov(:,:,timeIdx) = ExSq(:,:,timeIdx) - ...
                        xStats.mean(:,timeIdx)*xStats.mean(:,timeIdx)';

        % Correlation
        sigma =  diag(xStats.cov(:,:,timeIdx)).^0.5;
        sigmaInvTensor = 1./sigma * 1./sigma';
        xStats.corr(:,:,timeIdx) = xStats.cov(:,:,timeIdx) .* sigmaInvTensor;

    end

    alphaStats.var = EalphaSq - alphaStats.mean.^2;

    % Package raw data
    rawData.x = xStorage;
    rawData.u = uStorage;
    rawData.alpha = alphaStorage;
    rawData.rawStochasticValues = rawStochasticValues;
    rawData.Ygp = Ygp;
    rawData.Wcorrected = Wcorrected;
    rawData.simData = simDataStorage;
end