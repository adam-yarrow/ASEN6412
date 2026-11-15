function [Cmatrix, pceData] = FitPCE_ToMC(rawData_MC, p, times, fGP, d, ...
                        stochasticParam, stochasticParamGP, const)
    %{
        Fits a PCE to MC results.
        Solves a LSQ problem to fit psiMatrix*Cmatrix = x

        Assumptions:
            - Uniform random variables only
            - Using legendre polynomial basis functions

        Inputs:
            rawData_MC = structure with MC data for processing
            p = total order of the PCE basis    
            fGP = if true a gaussian process was also used, so need to process
            this separately
            d = dimension of random inputs
            stochasticParam = struct with stochastic input parameters
            stochasticParamGP = struct with GP parameters for CM_alpha
            const = model params struct

    %}

    %% Setup Vars
    nTimes = numel(times);
    N = numel(rawData_MC.simData); % Number of simulations available
    nStates = const.nStates;

    % Total basis functions
    P = factorial(p+d)/(factorial(p)*factorial(d)) - 1;

    % Storage Vessels
    Cmatrix = NaN(P+1, nTimes, nStates);
    RMSE_PCE = NaN(nStates, N);
    xEst = NaN(N, nTimes, nStates);
    psiMatrix = NaN(N, P+1);

    %% Get Input Variables
    % Standard Uniform Variables    
    [yNorm, ~] = GetNormalizedStdVariables(rawData_MC, stochasticParam);

    % Handle Gaussian Process Variables
    if fGP
        % assuming GP takes uses uniform random variables
        % assuming only a single GP present (CM_alpha)
        [yNormGP] = GetNormalizedGPVariables(rawData_MC, stochasticParamGP);

        % Make super Y
        yNorm = vertcat(yNorm, yNormGP);
    end
    
    %% Reshape State Matrix
    % Re-order raw states to be: N x nTimes x nStates
    xMatrix = permute(rawData_MC.x, [3,2,1]);

    %% Build psiMatrix once for all states
    [validBasisCombos] = GetValidPC_BasisOptions(d,p);
    for iMC = 1:N
        psiMatrix(iMC,:) = BuildPC_LegendreBasis(yNorm(:,iMC), validBasisCombos);
    end
    
    %% Solve System for each state across all times
    for iState = 1:nStates
        Cmatrix(:,:,iState) = psiMatrix \ xMatrix(:,:,iState);
        
        % Estimated value and RMSE calculation
        xEst(:,:,iState) =  psiMatrix*Cmatrix(:,:,iState);
        RMSE_PCE(iState,:) = sqrt(sum((xEst(:,:,iState) - xMatrix(:,:,iState)).^2, 2));
    end
    
    %% Plotting - Checking how good the solution is
    figure('Name','PCE RMSE Histograms');
    for iState = 1:nStates
        subplot(2,3,iState);
        histogram(RMSE_PCE(iState,:)*const.statePlottingSF(iState));
        grid on;
        xlabel(sprintf('%s (%s)',const.stateNames{iState}, const.statePlottingUnits{iState}))
        ylabel('Counts');
    end
    sgtitle('PCE RMSE Histograms');

    %% Save Relevant Data
    pceData = struct();
    pceData.psiMatrix = psiMatrix;
    pceData.xEst = permute(xEst, [3, 2, 1]); % Reshape to standard simData.x format
    pceData.RMSE = RMSE_PCE;
    pceData.basisOrders_jk = validBasisCombos; % This is the order of each univariate poly


end

%% helper functions
function [yNormGP] = GetNormalizedGPVariables(rawData_MC, stochasticParamGP)
    warning('this will probably need to change if i switch to uniform RV inputs');
    y = rawData_MC.rawY_GP; % dGP x N

    % Map y from [-K, K] -> [-1, 1] (from truncated gaussian process)
    K = stochasticParamGP.YtruncSigma;
    yNormGP = y / K; % Symmetric so can do this    
end

function [yNorm, dLocal] = GetNormalizedStdVariables(rawData_MC, stochasticParam)
    y = [rawData_MC.rawStochasticValues{:}]'; % d x N (N simulations)
    dLocal = size(y,1); % NOTE: only local cause GP handled separately
    
    % Map y from [a,b]-> [-1,1] for PCE basis
    yNorm = NaN(size(y));
    for iDim = 1:dLocal
        b = stochasticParam{iDim}.max;
        a = stochasticParam{iDim}.min;
    
        yNorm(iDim,:) =  2*(y(iDim,:) - (b+a)/2)/(b-a);
    end
end