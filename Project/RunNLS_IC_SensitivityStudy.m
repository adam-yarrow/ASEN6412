%{
    This function demonstrates the sensitivity to the initial NLS guess for
    CM_alpha coefficients.
%}


Startup();
%% Configuration
% Fix rng
rng(11112026,'twister');

%% Generate Trajectory Instance
% Real Simulation
mode_reality = 'linear-cmAlpha'; % Options = 'static-cmAlpha', 'linear-cmAlpha', 'Cmalpha-GP-staticMean'

% Running Model
[yMeas, t, nTimes, CM_alphaTrue, constReality, controller] = ...
    SetupInverseUQ(mode_reality);

%% NLSE Polynomial Config
% Bounds on CM_alpha_0 (e.g. the zeroth order term in a polynomial approx
% of CM_alpha)
CmAlphaMin = -2.5; % lower bound on search
CmAlphaMax = -1.5; % upper bound on search
nNLS_IC = 10;

% Estimator Settings
polyOrder = 1;
mode_estimator = 'poly-cmAlpha';

% Par-For
nWorkers = 8;

%% NLS IC Sensitivity
theta0_ConstRange = linspace(CmAlphaMin,CmAlphaMax, nNLS_IC);
nTerms = polyOrder + 1;

thetaEstStorage = NaN(nTerms, nNLS_IC);
PestStorage = NaN(nTerms, nTerms, nNLS_IC);
JnlsStorage = NaN(nNLS_IC);

tic;
WaitMessage = parfor_wait(nNLS_IC, 'Waitbar', true);
parfor (iThetaIC = 1:nNLS_IC, nWorkers)
    theta0 = [theta0_ConstRange(iThetaIC); zeros(nTerms-1,1)];
    [thetaEstStorage(:,iThetaIC), PestStorage(:,:,iThetaIC), nlsDetails] = ...
                    RunNLS(mode_estimator, theta0, ...
                        covInflationSF, yMeas, controller);

    JnlsStorage(iThetaIC) = nlsDetails.Jnls(end);
    WaitMessage.Send;
end
WaitMessage.Destroy
runTime = toc;

%% Plotting
figure('Name','Sensitivity of Jnls to IC');
plot(theta0_ConstRange, JnlsStorage);
grid on;
xlabel('CM_alpha Zeroth Order Parameter Initial Guess');
ylabel('J_{NLS}');
title('Sensitivity of Jnls to Initial Conditions');

% Plot diff realisations
tMatrix = BuildCmAlphaSensitivityMatrix(t, nTerms);

figure('Name','Inverse UQ - CM_Alpha Sensitivity to IC')
hold on;
plot(t, CM_alphaTrue,'k');
displayNames = cell(nNLS_IC + 1, 1);
displayNames{1} = 'True GP';
MSE_NLS_IC = NaN(nNLS_IC,1);
for i = 1:nNLS_IC
    CmAlpha_NLS = tMatrix * thetaEstStorage(:,i);
    MSE_NLS_IC(i) = sum((CmAlpha_NLS - CM_alphaTrue).^2);
    plot(t, CmAlpha_NLS,'--');
    displayNames{i+1} = sprintf('IC = %3.3f',theta0_ConstRange(i));
end
legend(displayNames);
grid on;
xlabel('Time (s)');
ylabel('Cm_{\alpha}');
title('Cm_{\alpha} vs Time for Different Initial Conditions');

% MSE vs IC
figure();
plot(theta0_ConstRange', MSE_NLS_IC);
grid on;
xlabel('CM_alpha Zeroth Order Parameter Initial Guess');
ylabel('MSE in CM_{\alpha}');
title('NLS MSE vs NLS Initial Guess Value');