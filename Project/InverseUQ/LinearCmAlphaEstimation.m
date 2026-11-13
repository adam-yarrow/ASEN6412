function [fTheta_Y, stats, rawData] = LinearCmAlphaEstimation(yMeas, rawTheta, ...
                                mode_estimator, mode_reality, ...
                                covInflationSF, fPlot, controller, nWorkers,...
                                integrationWeights, fPrior, thetaPaths)
    %{
        Estimates Cm_alpha given a set of measured data assuming a linear
        model.
    %}
                
    nTheta = size(rawTheta,1);

    % My Estimators Model
    constEstimator = BuildInverseUQParams(mode_estimator,covInflationSF);
         

    % Bayes me harder  
    [fTheta_Y, rawData] = DiscreteBayesFilter(rawTheta, thetaPaths, yMeas,...
                                                fPrior, ...
                                                controller, constEstimator,...
                                                nWorkers, ...
                                                integrationWeights);
    
    % Estimate Expected Values
    Etheta = sum(fTheta_Y .* rawTheta .* integrationWeights, 1)';
    
    EthetaSq = zeros(2,2);
    for i = 1:nTheta
        EthetaSq = EthetaSq + rawTheta(i,:)'*rawTheta(i,:) ...
            * fTheta_Y(i) * integrationWeights(i);
    end
    varTheta = EthetaSq - Etheta*Etheta';

    [~, MAPidx] = max(fTheta_Y);
    
    stats.MMSE = Etheta;
    stats.MAP = rawTheta(MAPidx,:)';
    stats.var = varTheta;
    
    % Plotting
    if fPlot
        figure('Name','Posterior of Cm_Alpha');
        scatter3(rawTheta(:,1), rawTheta(:,2), fTheta_Y);
        grid on;
        xlabel('Cm_{\alpha}');
        ylabel('d/dt(Cm_{\alpha})')
        zlabel('f_\theta(\theta)')
        title('Posterior distribution of Cm_{\alpha}', ...
            sprintf('Sim Model = %s - Estimator Model = %s',...
            mode_reality, mode_estimator));

    end
end