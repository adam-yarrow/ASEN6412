function [fY_Theta, fY_Theta_i, measInnov] = logLikelihoodFcn(y, theta, ...
                                            thetaPath, controller, const)
    %{  
        Calculates the likelihood of seeing a particular measurment given a
        particular dynamics model and measurement model.     

        Treats each measurement as independent. 

        Uses log likelihood to avoid issues with very small probabilities.

        Inputs:
            y = measurements: nMeas x nTimes
            theta = cell array of parameters to change for this likelihood
            calculation 
            thetaPath = cell array of paths to parameter locations in const
            controller = f(controller)
            const = baseline data structure.
    %}
    C = const.measModel.C;
    CxxEstimator = const.measModel.CxxEstimator;

    % Perturb constant structure
    constTemp = updateStruct(const, thetaPath, theta);

    % Run Model
    %% TODO - may replace this with a surrogate from fwd UQ if time/required
    simData = Simulation(controller, constTemp);

    % Meas Innov
    measInnov = y - C*simData.x;

    % Now every meas innov should be iid under error distribution
    [nMeas, nTimes] = size(y);

    fY_Theta = 0; % because log(1) = 0;
    n = const.nStates;
    fY_Theta_i = zeros(nTimes,1);
    for i = 1:nTimes
        % Noise is MV gaussian
        fY_Theta_i(i) = mvnpdf(measInnov(:,i), zeros(nMeas,1), CxxEstimator);

        % log(product of iid distributions) = sum(log(fX(xi))
        fY_Theta = fY_Theta + log(fY_Theta_i(i)); 
    end



end