function [stochasticParam] = ForwardUqPerturbationVars(mode)
    %{
        This function acts as a master source of truth for all forward UQ
        perturbations and their associated distributions.


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
    % Specify Random Variable Distributions
    CMalpha.path = 'aero.CM_alpha';
    CMalpha.distribution = 'trunc-gaussian';
    CMalpha.mean = 2.1335;
    CMalpha.sigma = 2.1335/10; % Assuming 10% variation per sigma
    CMalpha.truncatedSigmaLevel = 3; % truncate at +/- 3 sigma (+/- 30% variation)

    %% TODO - decide if changing these values too
    CMalphaSq.path = 'aero.CM_alphaSq';
    CMalphaSq.distribution = 'trunc-gaussian';
    CMalphaSq.mean = 6.2926;
    CMalphaSq.sigma = 6.2926/10; % Assuming 10% variation per sigma
    CMalphaSq.truncatedSigmaLevel = 3; % truncate at +/- 3 sigma (+/- 30% variation)

    Cm0.path = 'aero.CM_0';
    Cm0.distribution = 'trunc-gaussian';
    Cm0.mean = 1.8979E-1;
    Cm0.sigma = 1.8979E-2; % Assuming 10% variation per sigma
    Cm0.truncatedSigmaLevel = 3;
    
    %% IC Random Var Options
    V.path = 'xIC.V';
    V.distribution = 'trunc-gaussian'; %% TODO - change these to be trunc gaussian in future
    V.mean = 1500;
    V.sigma = 33; % 3-sigma = 100 m/s
    V.truncatedSigmaLevel = 3;
    
    theta.path = 'xIC.theta';
    theta.distribution = 'uniform'; % Proxy for alpha (if gamma = 0)
    theta.min = deg2rad(-5);
    theta.max = deg2rad(5);
    
    h.path = 'xIC.h';
    h.distribution = 'trunc-gaussian';
    h.mean = 28000;
    h.sigma = 500/3; % 3-sigma = 500 m/s
    h.truncatedSigmaLevel = 3;

    %% Process Different Options Based On Mode
    switch lower(mode)
        case lower('static-Cmalpha-only')
            stochasticParam = {CMalpha};

        case lower('static-Cmalpha-Cm0-CmalphaSq')
            stochasticParam = {CMalpha, CMalphaSq, Cm0};

        case lower('static-Cmalpha-wIC')
            stochasticParam = {CMalpha, V, theta, h};

        case lower('Cmalpha-GP-staticMean')
            % TODO
        case lower('Cmalpha-GP-linearMean')
            % TODO
        case lower('Cmalpha-GP-integralMean')
            % TODO

    end
end