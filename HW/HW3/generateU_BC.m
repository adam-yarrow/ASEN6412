function [uBC] = generateU_BC(Nmc, muTheta, sigmaTheta)
%{
    Generates boundary conditions for MC runs according to:
    u(0) = 0;
    u(1) ~ N(muTheta,sigmaTheta)

    uBC = Nmc x 2 vector of boundary conditions
%}
    uBC_1 = muTheta + sigmaTheta*randn(Nmc,1);
    uBC = [zeros(Nmc,1), uBC_1];
end