function [xStats, alphaStats, rawData] = StochasticColloc(const, Y, W, ...
                        controller, stochasticParam)
    %{
        This function runs stochastic collocation to estimate the
        trajectory statistics and forward UQ.

        const = baseline parameter set
        Y = d x N number of points to evaluate at
        W = N x 1 weights for each point
        controller = annonymous function representing the controller
        stochasticParam = cell array of structs defining how to perturb the
        parameters in the model.
    %}
    nAbsisscas = numel(W);

    % Storage
    nTimes = const.tEnd/const.dT + 1;
    xStorage = NaN(const.nStates, nTimes);
    uStorage = NaN(nTimes,1);
    alphaStorage = NaN(nTimes, 1);

    % Paths to perturb
    pathsToPerturb = cellfun(@(x) x.path, stochasticParam, 'UniformOutput', false);

    % Make Stochastic Colloc Inputs
    [stochasticValues, rawStochasticValues, Wcorrected] = MapSCCtoDistributions(...
                    stochasticParam, Y, W);
  
    % Run Stoch Colloc
    t0 = tic;
    wb = progressBar(1, nAbsisscas, t0, [], 'MC Runs');
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
            
        % Run Model
        simData = Simulation(controller, constTemp);
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
        progressBar(i, nAbsisscas, t0, wb);    
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
    rawData.simData = simDataStorage;
end