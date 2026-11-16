function [stochasticParam, stochasticParamGP, fGP] = ForwardUqPerturbationVars(mode, const)
    %{
        This function acts as a master source of truth for all forward UQ
        perturbations and their associated distributions.

        NOTE: if in gaussian process mode then this will return empty


        Example Types:
        Uniform
            data.path = 'dummy.struct.path';
            data.distribution = 'uniform';
            data.min = 0;
            data.max = 1;

        Truncated Gaussian
            data.path = 'dummy.struct.path';
            data.distribution = 'trunc-gaussian';
            data.mean = 0;
            data.sigma = 1;
            data.truncatedSigmaLevel = 3; % Symmetric truncation
    %}

    %% Aero Random Var Options
    CMalpha.path = 'aero.CM_alpha';
    CMalpha.distribution = 'uniform'; %'trunc-gaussian';
    CMalpha.min = -2.1335*1.3;
    CMalpha.max = -2.1335*0.7;
    
    %% IC Random Var Options
    V.path = 'xIC.V';
    V.distribution = 'uniform';%'trunc-gaussian'; 
    V.min = 1500 - 60;
    V.max = 1500 + 60;
    
    theta.path = 'xIC.theta';
    theta.distribution = 'uniform'; % Proxy for alpha (if gamma = 0)
    theta.min = deg2rad(-0.5);
    theta.max = deg2rad(0.5); % +/- 0.5 deg alpha pretty important
    
    h.path = 'xIC.h';
    h.distribution = 'uniform';
    h.min = 28000 - 100;
    h.max = 28000 + 100;

    %% Process Different Options Based On Mode
    fGP = false;
    stochasticParamGP = [];
    
    switch lower(mode)
        case lower('static-Cmalpha-only')
            stochasticParam = {CMalpha};

        case lower('static-Cmalpha-wIC')
            % stochasticParam = {CMalpha, V, theta, h};
            stochasticParam = {CMalpha, V, h, theta};

        case lower('Cmalpha-GP-staticMean')
            %{
                Constant mean and tuned parameters to try and make response
                look similar to static covariance case e.g. +/- 30% at
                3-sigma
            %}
            stochasticParamGP.meanFunction = @(t) -2.1335*ones(size(t)); % same as normal value
            targetSigma = 2.1335/10; % 10% variation in true
            sigmaSF = 1.0; % SF = 1 --> +-50% variation, SF = 0.75 --> +/- 30% total variation
            stochasticParamGP.sigma = targetSigma*sigmaSF; % tuned to give nice response (not too wide)
            stochasticParamGP.YtruncSigma = 3;
            stochasticParamGP.inputVariableMin = 0;
            stochasticParamGP.inputVariableMax = const.tEnd; % notinally time, but might change
            lSF = 0.32; % tuned to give approximately 10% MSE with d = 3;
            stochasticParamGP.l = const.tEnd*2*lSF; % Tuned to achieve around d = 3
            % stochasticParamGP.truncErrorMSE = 0.1; % 1% trunc error for KLE
            stochasticParamGP.d = 3;

            lambdaOptions = struct();
            lambdaOptions.mode = 'approximate';
            lambdaOptions.k = 200; % Assuming 200 eigenvalues beyond current eval adequate
            lambdaOptions.maxIts = 100;
            lambdaOptions.fPlot = false;

            stochasticParamGP.lambdaOptions = lambdaOptions;

            fGP = true;

            % Any static parameters can go here
            % TODO - for SOBOL indices add these parameters in as required
            stochasticParam = {V, h, theta};
    end
end