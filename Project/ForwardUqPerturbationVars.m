function [stochasticParam] = ForwardUqPerturbationVars()
    %{
        This function acts as a master source of truth for all forward UQ
        perturbations and their associated distributions.
    %}
    warning('How do i handle the fact distributions might be coupled?')

    % Specify Random Variable Distributions
    CMalpha.path = 'aero.CM_alpha';
    CMalpha.distribution = 'trunc-gaussian';
    CMalpha.mean = 2.1335;
    CMalpha.sigma = 2.1335/10; % Assuming 10% variation per sigma
    CMalpha.truncatedSigmaLevel = 3; % truncate at +/- 3 sigma (+/- 30% variation)

    %% TODO - decide if changing these values too
    CMalphaSq.path = 'aero.CM_alphaSq';
    CMalphaSq.distribution = 'uniform';
    CMalphaSq.mean = 6.2926;
    CMalphaSq.sigma = 6.2926/10; % Assuming 10% variation per sigma
    CMalphaSq.truncatedSigmaLevel = 3; % truncate at +/- 3 sigma (+/- 30% variation)

    Cm0.path = 'aero.CM_0';
    Cm0.distribution = 'uniform';
    Cm0.mean = 1.8979E-1;
    Cm0.sigma = 1.8979E-2; % Assuming 10% variation per sigma
    Cm0.truncatedSigmaLevel = 3;
    
    % V.distribution = 'uniform';
    % V.min = 1400;
    % V.max = 1600;
    V.path = 'xIC.V';
    V.distribution = 'trunc-gaussian'; %% TODO - change these to be trunc gaussian in future
    V.mean = 1500;
    V.sigma = 33; % 3-sigma = 100 m/s
    V.truncatedSigmaLevel = 3;
    
    theta.path = 'xIC.theta';
    theta.distribution = 'uniform'; % Proxy for alpha (if gamma = 0)
    theta.min = deg2rad(-5);
    theta.max = deg2rad(5);
    
    % h.distribution = 'uniform';
    % h.min = 27500;
    % h.max = 28500;
    h.path = 'xIC.h';
    h.distribution = 'trunc-gaussian';
    h.mean = 28000;
    h.sigma = 500/3; % 3-sigma = 500 m/s
    h.truncatedSigmaLevel = 3;

    % stochasticParam = {CMalpha, V, theta, h};
    
    % TODO - delete this once finished debugging
    stochasticParam = {CMalpha};
    warning('TODO - delete this temporary test fix');
end