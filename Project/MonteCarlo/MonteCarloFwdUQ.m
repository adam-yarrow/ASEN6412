function [xStats, alphaStats, rawData] = MonteCarloFwdUQ(const,Nmc,controller,...
    stochasticParam, stochasticParamGP)
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
    
    t0 = tic;
    wb = progressBar(1, Nmc, t0, [], 'MC Runs');
    simDataStorage = cell(Nmc, 1);
    Ygp = NaN(dGP, Nmc);

    for iMC = 1:Nmc
        % Overwrite perturbations - static variables
        constTemp = updateStruct(const, pathsToPerturb, stochasticValues{iMC});

        % Setup Gaussian process random variables
        if fGP
            Ktrunc = stochasticParamGP.YtruncSigma;
            Ygp(:,iMC) = trandn(-Ktrunc*ones(dGP,1), Ktrunc*ones(dGP,1));
            constTemp.aero.GP.Y = Ygp(:,iMC) ;
        end

        % Run Model
        simData = Simulation(controller, constTemp);
        simDataStorage{iMC} = simData;
        xStorage(:,:,iMC) = simData.x;
        uStorage(:,iMC) = simData.u;
        alphaStorage(:,iMC) = calcAlpha(simData.x);
    
        % Progress bar
        wb = progressBar(iMC, Nmc, t0, wb);
    end
    
    % Calculate Statistics
    [xStats, alphaStats] = CalculateMC_Statistics(xStorage, alphaStorage);
   
    % Package raw data
    rawData.x = xStorage;
    rawData.u = uStorage;
    rawData.alpha = alphaStorage;
    rawData.rawStochasticValues = rawStochasticValues;
    rawData.rawY_GP = Ygp;
    rawData.simData = simDataStorage;
end