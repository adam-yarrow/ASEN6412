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
    
    V.distribution = 'uniform'; %% TODO - change these to be trunc gaussian in future
    V.min = 1400;
    V.max = 1600;
    
    theta.distribution = 'uniform'; % Prox for alpha
    theta.min = deg2rad(-5);
    theta.max = deg2rad(5);
    
    h.distribution = 'uniform';
    h.min = 27500;
    h.max = 28500;
    
    stochasticParam = {CMq, V, theta, h};
    
end