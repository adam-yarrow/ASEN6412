%{
    Testing to see what i can get it to do
%}
Startup();

const = ModelParams();

xIC = [1500; 0; 0; 0; 28000]; % Guess
%% Basic Test
const.controller.sysIDActive = true;
PDControl = @(t,x,const) AlphaController(t,x,const, @(t) getReferenceAlpha(t));
simData = Simulation(xIC, PDControl, const);

const.aero.dynDeriv.model = 'static';
simData_CmqOn = Simulation(xIC, PDControl, const);

plotSimulationResults({simData,simData_CmqOn}, {'Cmq - OFF','Cmq - STATIC'});


%% Discrete Bayes filter for 1 variable on Basic Test
CmalphaMin = 0.8*const.aero.CM_alpha;
CmalphaMax = 1.2*const.aero.CM_alpha;
theta = linspace(CmalphaMin, CmalphaMax, 500); % uniform grid of theta

fTheta = ones(500,1) ./ (CmalphaMax - CmalphaMin); % Uniform prior

% Make measured data of single traj
Sxx = diag([10, deg2rad(5), deg2rad(5), deg2rad(5), 10]); % Dummy for now
CxxEstimator = 1*Sxx.^2; % Technically can be different
randVars = randn(size(simData.x));
y = simData.x + Sxx*randVars; % Added zero mean additive gaussian white noise

% Try Discrete Bayes Filter
[fTheta_Y, likelihoodVsTimeAndTheta, measInnov] = DiscreteBayesFilter(theta, y, fTheta, PDControl, tEnd, const, xIC, CxxEstimator);

figure();
plot(theta, fTheta_Y);
grid on;
xlabel('\theta');
ylabel('f_{\theta|Y}')


%% Bayesian estimation with random Cm_alpha
CmalphaMin = 0.5*const.aero.CM_alpha;
CmalphaMax = 1.5*const.aero.CM_alpha;
nTheta = 50;
theta = linspace(CmalphaMin, CmalphaMax, nTheta); % uniform grid of theta
Ntrials = 500;
constTemp = const;
sigmaCmAlpha = 0.2;
fThetaPrior = ones(nTheta,1) ./ (CmalphaMax - CmalphaMin);
fTheta_Yi = NaN(nTheta, Ntrials);
constTemp.tEnd = 5;
CxxEstimator = 1*Sxx.^2; % Allows for inflation if necessary
t0 = tic;
wb = progressBar(1, Ntrials, t0, [], 'MC Runs');
cmAlphaDistrib = const.aero.CM_alpha + sigmaCmAlpha*randn(Ntrials,1);
for i = 1:Ntrials    
    % Make Random CM Alpha
    constTemp.aero.CM_alpha = cmAlphaDistrib(i);
    simData_i = Simulation(xIC, PDControl, constTemp);

    % Make measurements
    randVars = randn(size(simData_i.x));
    y_i = simData_i.x + Sxx*randVars;

    % Bayes Filter
    fTheta_Yi(:,i) = DiscreteBayesFilter(theta, y_i, fThetaPrior, PDControl, tEnd, const, xIC, CxxEstimator);

    % Update prior with posterior
    % fThetaPrior = fTheta_Yi(:,i);

    wb = progressBar(i, Ntrials, t0, wb);
end

% Not sure if the best way to handle estimating the distribution of theta?
% Could try to recursively fit a distribution based on the individual
% estimates?

[~, MAPidx] = max(fTheta_Yi);
MAP = theta(MAPidx);

figure();
histogram(MAP,'Normalization','probability');


% Assuming Gaussian
muTheta = mean(MAP);
varTheta = var(MAP);


% Kind of like a gauss mixture?
averagePDF = sum(fTheta_Yi,2) ./  Ntrials;
%% IDEALLY Take this averaged pdf and use as a prior for the bayes filter above
%% Then loop this process until the averaged PDF converges

figure();
plot(theta, averagePDF);
hold on;
plot(theta, normpdf(theta, muTheta, varTheta.^0.5));
xlabel('Theta');
ylabel('f_theta');
legend('Averaged PDF','Gaussian MAP Approximation')



% figure();
% plot(theta, fTheta_Yi(:,end));
% grid on;
% xlabel('\theta');
% ylabel('f_{\theta|Y}')


%% Testing Stochastic Param
PDControl = @(t,x,const) AlphaController(t,x,const);

const = ModelParams();


Nsims = 1000;

tEnd = 10;
Ndata = tEnd/const.dT + 1;
times = 0:const.dT:tEnd;

% Random Variables
% TODO - maybe replace with beta distributions?
sigmaCmAlpha = 0.1;
% sigmaCmAlphaSq = 0.6;
% sigmaCmDeltaE = 0.1;
% sigmaCM0 = 0.005;

% CM0 = const.aero.CM_0 + sigmaCM0*randn(Nsims,1);
CmAlpha = const.aero.CM_alpha + sigmaCmAlpha*randn(Nsims, 1);
% CmAlphaSq = const.aero.CM_alphaSq + sigmaCmAlphaSq*randn(Nsims, 1);
% CmDeltaE = const.aero.CM_deltaE  + sigmaCmDeltaE*randn(Nsims,1);

gamma0 = 0;
muAlpha0 = 0;
sigmaAlpha0 = deg2rad(10/3);
alpha0 = muAlpha0 + sigmaAlpha0*randn(Nsims,1);

constTemp = const;
xICTemp = xIC;
 
% Running Sims
deltaE_N = NaN(Ndata, Nsims);
x_N = NaN(const.nStates, Ndata, Nsims);
alpha_N = NaN(Ndata, Nsims);


for i = 1:Nsims
    % Make Random Vars
    constTemp.aero.CM_alpha = CmAlpha(i);
    % constTemp.aero.CM_alphaSq = CmAlphaSq(i);
    % constTemp.aero.CM_deltaE = CmDeltaE(i);
    % constTemp.aero.CM_0 = CM0(i);

    xICTemp(4) = alpha0(i) + gamma0;

    simDataTemp = Simulation(xICTemp, PDControl, constTemp);
    
    % Pull Useful data
    deltaE_N(:,i) = simDataTemp.u(2,:);
    x_N(:,:,i) = simDataTemp.x;
    alpha_N(:,i) = calcAlpha(simDataTemp.x);
end

deltaE_Mean = mean(deltaE_N,2);
deltaE_Var = var(deltaE_N,[],2);

figure();
plot(times, rad2deg(alpha_N));

figure(); 
hold on; 
plot(times, deltaE_Mean); 
plot(times, deltaE_Mean + sqrt(deltaE_Var)); 
plot(times, deltaE_Mean - sqrt(deltaE_Var));





%% Discrete Bayes Filter Functions
% now need likelihood function
function [fY_Theta, fY_Theta_i, measInnov] = fYgivenTheta(y, theta_i, PDControl, tEnd, const, xIC, CxxEstimator)
    % Log likelihood
  
    % Run simulation using theta_i parameter
    const.aero.CM_alpha = theta_i; 
    [simData] = Simulation(xIC, PDControl, const); % Assuming you know xIC
    measInnov = y - simData.x;

    % Now every meas innov should be iid under error distribution
    NTimes = size(y,2);
    fY_Theta = 0;
    n = const.nStates;
    fY_Theta_i = zeros(NTimes,1);
    for i = 1:NTimes
        % Noise is MV gaussian
        fY_Theta_i(i) = 1./((2*pi)^(n/2)*sqrt(det(CxxEstimator)))  * ...
                exp(-1/2*(measInnov(:,i)'*inv(CxxEstimator)*measInnov(:,i)));
        fY_Theta = fY_Theta + log(fY_Theta_i(i)); % product of iid distributions
    end
       
end

% Discrete Bayes Filter for a single trajectory
function [fTheta_Y, likelihoodVsTimeAndTheta, measInnov] = DiscreteBayesFilter(theta, y, fThetaPrior, PDControl, tEnd, const, xIC, CxxEstimator)
    nTheta = numel(theta);
    fTheta_Y = zeros(nTheta,1);
    likelihoodVsTimeAndTheta = NaN(size(y,2), nTheta);
    measInnov = NaN(size(y,1), size(y,2), nTheta);
    for i = 1:nTheta
        [logLikelihood, likelihoodVsTimeAndTheta(:,i), measInnov(:,:,i)] = fYgivenTheta(y, theta(i), PDControl, tEnd, const, xIC, CxxEstimator);
        fTheta_Y(i) = logLikelihood + log(fThetaPrior(i));
    end

    % Convert From Log to Normal Domain
    fTot = LogSumExp(fTheta_Y);
    fTheta_Y = exp(fTheta_Y - fTot);

    % Normalize into a Proper Density (e.g. sums to 1.0)
    fTheta_Y = fTheta_Y ./ trapz(theta, fTheta_Y);

    % % Old normalization
    % fTheta_Y = fTheta_Y - max(fTheta_Y); % arbitrarily normalize to make exp less painful
    % fY = sum(exp(fTheta_Y));
    % 
    % fTheta_Y = exp(fTheta_Y) ./ fY; 
end