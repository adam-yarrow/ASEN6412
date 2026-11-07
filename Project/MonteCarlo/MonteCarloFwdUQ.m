function [xStats, alphaStats, rawData] = MonteCarloFwdUQ(const,Nmc,controller,...
    pathsToPerturb,stochasticParam)
    %{
        This function runs a MC simulation to estimate the trajectory
        statistics and forward UQ.
        
        const = baseline parameter set
        Nmc = Number of MC simulations to run
        controller = annonymous function representing the controller
        pathsToPerturb = cell array containing the paths to the variables
        in const to perturb stochastically
        stochasticParam = cell array of structs defining how to perturb the
        parameters in the model.
    %}
    
    % Storage
    nTimes = const.tEnd/const.dT + 1;
    xStorage = NaN(const.nStates, nTimes, Nmc);
    uStorage = NaN(nTimes, Nmc);
    alphaStorage = NaN(nTimes, Nmc);
    
    % Generate Random Inputs for each MC case
    [stochasticValues, rawStochasticValues] = DrawRandomVariables(stochasticParam, Nmc);
    
    t0 = tic;
    wb = progressBar(1, Nmc, t0, [], 'MC Runs');
    simDataStorage = cell(Nmc, 1);
    for iMC = 1:Nmc
        % Overwrite perturbations
        constTemp = updateStruct(const, pathsToPerturb, stochasticValues{iMC});
    
        % Run Model
        simData = Simulation(controller, constTemp);
        simDataStorage{iMC} = simData;
        xStorage(:,:,iMC) = simData.x;
        uStorage(:,iMC) = simData.u;
        alphaStorage(:,iMC) = calcAlpha(simData.x);
    
        % Progress bar
        progressBar(iMC, Nmc, t0, wb);
    end
    
    % Calculate Statistics
    [xStats, alphaStats] = CalculateMC_Statistics(xStorage, alphaStorage);
   
    % Package raw data
    rawData.x = xStorage;
    rawData.u = uStorage;
    rawData.alpha = alphaStorage;
    rawData.rawStochasticValues = rawStochasticValues;
    rawData.simData = simDataStorage;
end