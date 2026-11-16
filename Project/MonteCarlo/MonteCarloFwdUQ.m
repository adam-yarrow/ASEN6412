function [xStats, alphaStats, rawData] = MonteCarloFwdUQ(const,Nmc,controller,...
    stochasticParam, stochasticParamGP, nWorkers)
    %{
        This function runs a MC simulation to estimate the trajectory
        statistics and forward UQ.
        
        const = baseline parameter set
        Nmc = Number of MC simulations to run
        controller = annonymous function representing the controller
        stochasticParam = cell array of structs defining how to perturb the
        parameters in the model.
        stochasticParamGP = data for stochastic gaussian process of
        Cm_alpha
    %}
    
    % Storage
    nTimes = const.tEnd/const.dT + 1;
    xStorage = NaN(const.nStates, nTimes, Nmc);
    uStorage = NaN(nTimes, Nmc);
    alphaStorage = NaN(nTimes, Nmc);

    % Paths to perturb
    pathsToPerturb = cellfun(@(x) x.path, stochasticParam, 'UniformOutput', false);
    
    % Generate Random Inputs for each MC case
    [stochasticValues, rawStochasticValues] = DrawRandomVariables(stochasticParam, Nmc);

    % Setup Gaussian Process Cm_alpha
    fGP = false;
    dGP = 0;
    if ~isempty(stochasticParamGP)
        const.aero.GP.fActive = true;
        const.aero.GP.dataKLE = stochasticParamGP.dataKLE;
        dGP = stochasticParamGP.dataKLE.d;
        fGP = true;
    end
    
    simDataStorage = cell(Nmc, 1);
    Ygp = NaN(dGP, Nmc);

    WaitMessage = parfor_wait(Nmc, 'Waitbar', true);
    parfor (iMC = 1:Nmc, nWorkers)
        % Overwrite perturbations - static variables
        constTemp = updateStruct(const, pathsToPerturb, stochasticValues{iMC});

        % Setup Gaussian process random variables. Assuming underlying
        % uniform distribution for Ygp = [-Ktrunc, Ktrunc]
        if fGP
            Ktrunc = stochasticParamGP.YtruncSigma;
            Ygp(:,iMC) = Ktrunc*rand(dGP,1)
            constTemp.aero.GP.Y = Ygp(:,iMC) ;
        end

        % Run Model
        simData = Simulation(controller, constTemp);
        simDataStorage{iMC} = simData;
        xStorage(:,:,iMC) = simData.x;
        uStorage(:,iMC) = simData.u;
        alphaStorage(:,iMC) = calcAlpha(simData.x);
    
        % Progress bar
        WaitMessage.Send;
    end

    WaitMessage.Destroy
    
    % Calculate Statistics
    % [xStats, alphaStats] = CalculateMC_Statistics(xStorage, alphaStorage);
    [xStats, alphaStats] = CalculateMC_StatisticsFast(xStorage, alphaStorage); % 50x's faster than my implementation
 
   
    % Package raw data
    rawData.x = xStorage;
    rawData.u = uStorage;
    rawData.alpha = alphaStorage;
    rawData.rawStochasticValues = rawStochasticValues;
    rawData.rawY_GP = Ygp;
    rawData.simData = simDataStorage;
end