function [yNormGP] = GetNormalizedGPVariables(rawData_MC, stochasticParamGP)
    y = rawData_MC.Ygp; % dGP x N

    % Map y from [-K, K] -> [-1, 1] (from truncated gaussian process)
    K = stochasticParamGP.YtruncSigma;
    yNormGP = y / K; % Symmetric so can do this    
end
