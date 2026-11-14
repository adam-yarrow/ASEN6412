function [theta_NLS, P_NLS, nlsDetails] = RunNLS(mode_estimator, theta0, ...
                                        covInflationSF, yMeas, controller,...
                                        fPrintUpdates)
    if nargin < 6
        fPrintUpdates = true;
    end
    % Pull out options
    nlsOptions = ModelParams('nlsOptions');
    
    % Setup Solver Inputs
    constEstimator = BuildInverseUQParams(mode_estimator,covInflationSF);
    HFunc = @(thetaK) CreateHmatrix(thetaK, constEstimator, controller);
    hFunc = @(thetaK) CreateMeasEstimates(thetaK, constEstimator, controller);
    
    % Run Solver
    nTimes = size(yMeas,2);
    R_NLS = diag(repmat(diag(constEstimator.measModel.Cxx)*covInflationSF, nTimes, 1));
    yMeasNLS = reshape(yMeas, [], 1);
    [theta_NLS, P_NLS, nlsDetails] = NLS(theta0,hFunc,HFunc,...
                yMeasNLS,R_NLS,nlsOptions,fPrintUpdates);
end