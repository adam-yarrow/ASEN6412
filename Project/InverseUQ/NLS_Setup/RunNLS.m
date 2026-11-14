function [theta_NLS, P_NLS, nlsDetails] = RunNLS(mode_estimator, theta0, ...
                                        covInflationSF, yMeas, controller,...
                                        fPrintUpdates)
    if nargin < 6
        fPrintUpdates = true;
    end
    % Options Struct for NLS Solver
    nlsOptions = struct();
    nlsOptions.maxIterations = 15;
    nlsOptions.JnlsRelTol =  1E-12;
    nlsOptions.x0Tol = 1E-8;
    nlsOptions.maxAlphaIterations = 50;
    nlsOptions.alphaRelTol = 1E-4; % Relative change between estimates of optimal alpha that triggers convergence
    nlsOptions.minAlpha = 1E-10;
    nlsOptions.fAlphaGoldenSearch = true;
    nlsOptions.alphaSF = 0.5;
    
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