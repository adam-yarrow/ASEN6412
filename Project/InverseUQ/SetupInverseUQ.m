function [yMeas, t, nTimes, CM_alphaTrue, constReality, controller] = ...
                SetupInverseUQ(mode_reality)
    %{
        Builds variables common to all inverse UQ scripts.
    %}
    
    % Model Baseline
    controller =  @(t,x,const) AlphaController(t, x, const, @(t) getReferenceAlpha(t));

    constReality = BuildInverseUQParams(mode_reality);
    
    [simData] = Simulation(controller, constReality);
    yMeas = getMeasurements(constReality, simData.x);
    nTimes = size(yMeas,2);
    t = simData.times';
    
    simDataNoisy = simData;
    simDataNoisy.x = yMeas;
    
    % Plot truth and measurements
    plotSimulationResults({simData, simDataNoisy},{'Truth Trajectory','Measurements'});
    
    switch lower(mode_reality)
        case lower('static-cmAlpha')
            CM_alphaTrue = ones(nTimes,1)*constReality.aero.CM_alpha;

        case lower('linear-cmAlpha')
            CM_alphaTrue = constReality.aero.CM_alpha + ...
                t*constReality.aero.movMeanCM_alpha.deriv;
        
        case lower('Cmalpha-GP-staticMean')
            phi = MakePhiMatrix(t, constReality.aero.GP.dataKLE);   
            CM_alphaTrue = getRV_KLE_Fast(t, constReality.aero.GP.dataKLE, ...
                                    phi, constReality.aero.GP.Y);
    end
end