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
    stochasticValues = {};
    rawStochasticValues = {};


    [stochasticValues, rawParams, Wcorrected] = MapSCCtoDistributions(...
                    stochasticParams, Y, W)
    %TODO - make this


    % Run Stoch Colloc
    t0 = tic;
    wb = progressBar(1, Nmc, t0, [], 'MC Runs');
    simDataStorage = cell(nAbsisscas, 1);

    for i = 1:nAbsisscas
        % Overwrite perturbations
        constTemp = updateStruct(const, pathsToPerturb, stochasticValues{iMC});
            
        % Run Model
        simData = Simulation(controller, constTemp);
        simDataStorage{iMC} = simData;
        xStorage(:,:,iMC) = simData.x;
        uStorage(:,iMC) = simData.u;
        alphaStorage(:,iMC) = calcAlpha(simData.x);

        % Calculate Statistics Using Stochastic Collocation
        
        %% TODO - work this out.

    
        % Progress bar
        progressBar(iMC, Nmc, t0, wb);
    
    end



    % Package raw data
    rawData.x = xStorage;
    rawData.u = uStorage;
    rawData.alpha = alphaStorage;
    rawData.rawStochasticValues = rawStochasticValues;
    rawData.simData = simDataStorage;


end