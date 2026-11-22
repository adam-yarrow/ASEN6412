function [fTheta_Y, stats, rawData] = CmAlphaEstimation(yMeas, rawTheta, ...
                                mode_estimator, mode_reality, ...
                                covInflationSF, fPlot, controller, nWorkers)
    %{
        Estimates Cm_alpha given a set of measured data.
    
        TODO - describe more
    %}
    
    nTheta = size(rawTheta,1);

    % Generate Trapezoidal Integration weights - using trapz as have discrete
    % grid that is fairly sparse
    dTheta = rawTheta(2) - rawTheta(1);
    integrationWeightsTheta = dTheta/2 * [1; 2*ones(nTheta-2,1); 1]; % TrapZ weights
    
    % Prior - Uniform
    thetaMax = max(rawTheta);
    thetaMin = min(rawTheta);
    fThetaPrior = ones(size(rawTheta))*(1/(thetaMax - thetaMin));
    
    % My Estimators Model
    constEstimator = BuildInverseUQParams(mode_estimator,covInflationSF);
    
    switch lower(mode_estimator)
        case lower('static-cmAlpha')
            thetaPath = {'aero.CM_alpha'};  
            theta = rawTheta;
        case lower('linear-cmAlpha')
            thetaPath = {'aero.CM_alpha'};  
            theta = rawTheta;
        case lower('Cmalpha-GP-staticMean')
            thetaPath = {'aero.GP.meanX'};
            theta = cell(nTheta,1);
            for iTheta = 1:nTheta
                theta{iTheta} = @(t) rawTheta(iTheta)*ones(size(t));
            end
        case lower('Cmalpha-GP-movMean')
            thetaPath = {'aero.GP.meanX'};
    
            theta = cell(nTheta,1);
            for iTheta = 1:nTheta
                theta{iTheta} = @(t) rawTheta(iTheta)*ones(size(t));
            end
    end
    
    [fTheta_Y, rawData] = DiscreteBayesFilter(theta, thetaPath, yMeas,...
                                                fThetaPrior, ...
                                                controller, constEstimator,...
                                                nWorkers, ...
                                                integrationWeightsTheta);
    
    % Estimate Expected Values
    Etheta = sum((fTheta_Y .* rawTheta) .* integrationWeightsTheta);
    EthetaSq = sum((fTheta_Y .* rawTheta .^2) .* integrationWeightsTheta);
    varTheta = EthetaSq - Etheta^2;

    [~, MAPidx] = max(fTheta_Y);
    
    stats.MMSE = Etheta;
    stats.MAP = rawTheta(MAPidx);
    stats.var = varTheta;
    
    % Plotting
    if fPlot
        figure('Name',sprintf('Posterior of Cm_Alpha - Sim Model = %s - Estimator Model = %s',...
            mode_reality, mode_estimator));
        plot(rawTheta, fTheta_Y);
        grid on;
        xlabel('Cm_{\alpha}');
        ylabel('f_\theta(\theta)')
        title('Posterior distribution of Cm_{\alpha}', ...
            sprintf('Sim Model = %s - Estimator Model = %s',...
            mode_reality, mode_estimator));

    end
end