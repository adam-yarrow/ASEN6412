function [fTheta_Y, stats, rawData] = LinearCmAlphaEstimation(yMeas, rawTheta, ...
                                mode_estimator, mode_reality, ...
                                covInflationSF, fPlot, controller, nWorkers,...
                                integrationWeights, fPrior)
    %{
        Estimates Cm_alpha given a set of measured data assuming a linear
        model.
    %}
    
    nTheta = size(rawTheta,1);
        
    
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
        figure('Name','Posterior of Cm_Alpha');
        plot(rawTheta, fTheta_Y);
        grid on;
        xlabel('Cm_{\alpha}');
        ylabel('f_\theta(\theta)')
        title('Posterior distribution of Cm_{\alpha}', ...
            sprintf('Sim Model = %s - Estimator Model = %s',...
            mode_reality, mode_estimator));

    end
end