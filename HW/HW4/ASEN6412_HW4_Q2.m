%{
    TODO - readme
%}
%% Startup
addpath('../KLE');

% Fix random seed for repeatability
rng(30102026);

%% Fixed Parameters
% G(x,w) 
muK = 1.0;
sigmaKxx = 0.1;
lKxx = 2;

% Truncation Terms
d = 2;

% Domain
xMax = 1;
xMin = 0;
Nx = 101; % From Previous HW3 solution
x = linspace(xMin, xMax, Nx)'; % Fixed domain

% KLE Random Var Uniform Domain
Ya = [-1, -1];
Yb = [1, 1];

% Boundary Conditions
uBC = [0, 0]; % Corresponds to x = [0, 1]


%% Setup KLE Expansion Common Parameters
% Get d eigenfunctions
lambdaKxx = NaN(d,1);
omegaKxx = NaN(d,1);
aKxx = (xMax - xMin)/2; % Calculate domain boundary

for i = 1:d
    omegaKxx(i) = getOmegaKLE(i,aKxx,lKxx);
    lambdaKxx(i) = getLambdaKLE(lKxx,omegaKxx(i),sigmaKxx);
end

% Generate KLE Data
meanKFunc = @(x) muK * ones(size(x));
dataKLE = generateKLE_Struct(lambdaKxx, omegaKxx, meanKFunc, aKxx, d);

% Generate Phi Matrix to save time for our grid
dx = x(2) - x(1);
x_itHalf = x - dx/2; 
x_itHalf = x_itHalf(2:end); % First entry outside of domain -> remove

% Pre-make eigenfunctions - I AM SPEED!
% Eigenfunctions don't change between grids
[phiMatrix] = MakePhiMatrix(x_itHalf, dataKLE);

dataKLE.x_itHalf = x_itHalf;
dataKLE.phiMatrix = phiMatrix;

%% TODO - be careful because MonteCarloPDE assumes Gxx = exp(Kxx)
warning('fix me');


%% LHS Solution - <u(x)> and var(u(x))
N_LHS = [1E1, 1E2, 1E3, 1E4, 1E5];
NumLHS = numel(N_LHS);

fPlotLHS = false;

for iLHS = 1:NumLHS
    Y_LHS = getLHS_SamplesUniformDist(N_LHS(iLHS), Ya, Yb);     
    [u, uStats] = MonteCarloPDE_HW4(x, uBC, N_LHS(iLHS), dataKLE, Y_LHS, fPlotLHS);

    warning('TOdo - save data somewhere, and actually calculate stats')
end
