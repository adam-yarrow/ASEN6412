function [pathsToPerturb, stochasticParam] = ForwardUqPerturbationVars()
    %{
        This function acts as a master source of truth for all forward UQ
        perturbations and their associated distributions.
    %}
    % Everything must appear in the same order it was built in the pathsToPerturb
    % variable
    pathsToPerturb = {'aero.dynDeriv.CM_q_static',...
                       'xIC.V',...
                       'xIC.theta',...
                       'xIC.h'};

    warning('How do i handle the fact distributions might be coupled?')

    % Specify Random Variable Distributions
    CMq.distribution = 'uniform';
    CMq.min = -0.2*100;
    CMq.max = -0.2*500;
    
    % V.distribution = 'uniform';
    % V.min = 1400;
    % V.max = 1600;
    V.distribution = 'trunc-gaussian'; %% TODO - change these to be trunc gaussian in future
    V.mean = 1500;
    V.sigma = 33; % 3-sigma = 100 m/s
    V.truncatedSigmaLevel = 3;
    
    theta.distribution = 'uniform'; % Proxy for alpha (if gamma = 0)
    theta.min = deg2rad(-5);
    theta.max = deg2rad(5);
    
    % h.distribution = 'uniform';
    % h.min = 27500;
    % h.max = 28500;
    h.distribution = 'trunc-gaussian';
    h.mean = 28000;
    h.sigma = 500/3; % 3-sigma = 500 m/s
    h.truncatedSigmaLevel = 3;

    stochasticParam = {CMq, V, theta, h};
    
end