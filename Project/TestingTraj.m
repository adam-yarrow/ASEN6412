%{
    Testing to see what i can get it to do
%}

const = ModelParams();

tEnd = 50;
xIC = [1500; 0; 0; 0; 28000]; % Guess
%% Basic Test
PDControl = @(t,x,const) BasicController(t,x,const);
simData = Simulation(xIC, PDControl, tEnd, const);
plotSimulationResults({simData}, {'test'});


%% Testing Stochastic Param
PDControl = @(t,x,const) BasicController(t,x,const);

const = ModelParams();


Nsims = 1000;

tEnd = 10;
Ndata = tEnd/const.dT + 1;
times = 0:const.dT:tEnd;

% Random Variables
% TODO - maybe replace with beta distributions?
sigmaCmAlpha = 0.1;
sigmaCmAlphaSq = 0.6;
sigmaCmDeltaE = 0.1;
sigmaCM0 = 0.005;

CM0 = const.aero.CM_0 + sigmaCM0*randn(Nsims,1);
CmAlpha = const.aero.CM_alpha + sigmaCmAlpha*randn(Nsims, 1);
CmAlphaSq = const.aero.CM_alphaSq + sigmaCmAlphaSq*randn(Nsims, 1);
CmDeltaE = const.aero.CM_deltaE  + sigmaCmDeltaE*randn(Nsims,1);

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
    constTemp.aero.CM_alphaSq = CmAlphaSq(i);
    constTemp.aero.CM_deltaE = CmDeltaE(i);
    constTemp.aero.CM_0 = CM0(i);

    xICTemp(4) = alpha0(i) + gamma0;

    simDataTemp = Simulation(xICTemp, PDControl, tEnd, constTemp);
    
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
plot(times, deltaE_Mean - sqrt(deltaE_Var))

%% Controller
function [u, exit] = BasicController(t,x,const)    
    exit = 0;

    f = 1/10;
    omega = 2*pi*f;
    A = deg2rad(10);
    alphaDes = A * sin(omega*t);
    alphaRateDes = -A*omega*cos(omega*t);

    % Do alpha control
    Kd = 20;
    Kp = 60; % just a guess lol

    alpha = calcAlpha(x);

    % PD Controller
    virtualElevator =  (-Kd*(x(3) - alphaRateDes) - Kp*(alpha - alphaDes));

    % Saturation
    deltaE_cmd = min(max(-virtualElevator, const.constraint.deltaEBounds(1)), ...
                    const.constraint.deltaEBounds(2)); 
        
    % Build Inputs
    u = [0; deltaE_cmd; 0]; % No engine and no canard
end