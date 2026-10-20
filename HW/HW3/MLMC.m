%{
    MLMC bonus question.

    This script shows the relative cost of MLMC vs MC methods for different
    MSE levels and grid refinement.
    
%}
StartupHW3();

%% Fixed Parameters
% G(x,w) 
muG = 1.0;
sigmaGxx = 2;
lGxx = 0.2;

% Truncation Terms
d = 10;

% x = 1 BC: theta(w) ~ N(mu,sigma)
muTheta = 0;
sigmaTheta = 0.1;

xMax = 1;
xMin = 0;

%% Setup KLE Expansion Common Parameters
% Get d eigenfunctions
lambdaGxx = NaN(d,1);
omegaGxx = NaN(d,1);
aGxx = (xMax - xMin)/2; % Calculate domain boundary

for i = 1:d
    omegaGxx(i) = getOmegaKLE(i,aGxx,lGxx);
    lambdaGxx(i) = getLambdaKLE(lGxx,omegaGxx(i),sigmaGxx);
end

% Generate KLE Data
meanGFunc = @(x) muG * ones(size(x));
dataKLE = generateKLE_Struct(lambdaGxx, omegaGxx, meanGFunc, aGxx, d);

fPlotMC = false;

%% Generate Pilot Sample Estimates
Npilot = 100;

NxGrid = [8, 32, 128, 1024];
l = 0:3;
nGrid = numel(NxGrid);

VarUestimate_Ql_xHalf = NaN(nGrid,1);
VarUestimates_Qdiff_xHalf = NaN(nGrid-1,1); % [Var(Q1-Q0), ...]

% Generate new independent random raws on BC and KLE Y(w) for each grid
% size
for i = 1:nGrid   
    if i == 1 % level 0 case
        Y = randn(d, Npilot);
        [uBC] = generateU_BC(Npilot, muTheta, sigmaTheta);

        x = linspace(xMin,xMax,NxGrid(i))';
        [~, uStats] = MonteCarloPDE(x, uBC, dataKLE, Y, fPlotMC);
        
        VarUestimate_Ql_xHalf(i) = interp1(x, uStats.var(:,end), 0.5);


    else % other levels case
        
        Y = randn(d, Npilot);
        [uBC] = generateU_BC(Npilot, muTheta, sigmaTheta);

        x_l = linspace(xMin,xMax,NxGrid(i))';
        [u_l, uStats_l] = MonteCarloPDE(x_l, uBC, dataKLE, Y, fPlotMC);
        VarUestimate_Ql_xHalf(i) = interp1(x_l, uStats_l.var(:,end), 0.5); % For MC Calcs

        x_lt1 = linspace(xMin,xMax,NxGrid(i-1))';
        [u_lt1, uStats_lt1] = MonteCarloPDE(x_lt1, uBC, dataKLE, Y, fPlotMC);

        Qdiff = interp1(x_l, u_l, 0.5) - interp1(x_lt1, u_lt1, 0.5);
        VarUestimates_Qdiff_xHalf(i-1) = var(Qdiff);
    end
end

%% Now find optimal N and cost (C)
epsilon = [0.1, 0.01, 0.001];
nEps = numel(epsilon);

relativeCostPerGrid = (NxGrid)' / (NxGrid(1)); % E.g. Normalized cost to smallest grid size

Ctot_MC = NaN(nEps, nGrid);
Ntot_MC = NaN(nEps, nGrid);

Ctot_MLMC = NaN(nEps, 1);

Nl_MLMC = NaN(nEps, nGrid);

for k = 1:nEps
    epsK = epsilon(k);
   [Ctot_MC(k,:), Ntot_MC(k,:)] = getMC_Cost(relativeCostPerGrid, ...
                            VarUestimate_Ql_xHalf, epsK);

    [Ctot_MLMC(k), ...
        Nl_MLMC(k,1), ...
        Nl_MLMC(k,2:end)] = getMLMC_Cost(relativeCostPerGrid, ...
                            VarUestimate_Ql_xHalf, VarUestimates_Qdiff_xHalf,...
                            epsK);
    
end

%% Plotting
figure('Name','MC vs MLMC Cost');
loglog(epsilon, Ctot_MLMC,'k--','DisplayName','MLMC','LineWidth',1.5);
hold on;
for i = 1:nGrid
    loglog(epsilon, Ctot_MC(:,i),'DisplayName',sprintf('MC Level = %i',l(i)),...
        'LineWidth',1.5);
end
grid on;
legend();
ylabel("Total Cost");
xlabel('\epsilon')
title('Cost Comparison of MLMC and MC Methods')


%% Helper Functions
function [Ctot, Ntot] = getMC_Cost(relativeCostPerGrid, Var_Ql, epsilon_k)
    Ntot = ceil(2*Var_Ql / epsilon_k^2);
    Ctot = Ntot .* relativeCostPerGrid;
end

function [Ctot, N0min, Nlmin] = getMLMC_Cost(relativeCostPerGrid, Var_Ql, Var_Qdiff, epsilon_k)
    VarQ0 = Var_Ql(1);
    C0 = relativeCostPerGrid(1);
    Cl = relativeCostPerGrid(2:end) + relativeCostPerGrid(1:end-1); % E.g. because running two diff size grids
    C = sqrt(VarQ0*C0) + sum(sqrt(Var_Qdiff .* Cl));

    N0min = 2*C/epsilon_k^2 * sqrt(VarQ0/C0);
    Nlmin = 2*C/epsilon_k^2 * sqrt(Var_Qdiff ./ Cl);
    
    % Round up to nearest integer multiples of N (such that achieve >=
    % condition)
    Nlmin = round(Nlmin);
    N0min = round(N0min);

    Ctot = N0min*C0 + dot(Cl, Nlmin);    
end
