%{
    TODO - readme
%}
%% Startup
addpath(genpath('..\..\Library\'));

% Fix random seed for repeatability
rng(01112026);

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


%% LHS Solution - <u(x)> and var(u(x))
N_LHS = [1E1, 1E2, 1E3, 1E4, 1E5];
NumLHS = numel(N_LHS);

fPlotLHS = false;

uLHS = cell(NumLHS,1);
uStatsLHS = cell(NumLHS,1);

mu_LHS = NaN(Nx, NumLHS);
var_LHS = NaN(Nx, NumLHS);

for iLHS = 1:NumLHS
    Y_LHS = getLHS_SamplesUniformDist(N_LHS(iLHS), Ya, Yb);     
    [uLHS{iLHS}, uStatsLHS{iLHS}] = MonteCarloPDE_HW4(x, uBC, N_LHS(iLHS), ...
                            dataKLE, Y_LHS);
    mu_LHS(:,iLHS) = uStatsLHS{iLHS}.mu;
    var_LHS(:,iLHS) = uStatsLHS{iLHS}.var;

end

% Plotting
plotUxVsNpoints('LHS', x, mu_LHS, var_LHS, N_LHS, true);
plotUxVsNpointsDifference('LHS', x, mu_LHS, var_LHS, N_LHS, true);

%% Stochastic Collocation - with CC Rule
% Reference Solution
levels = [0, 1, 2, 3, 4, 5];
lRefIdx = 6;
N_levels = [1, 2.^levels(2:end) + 1];
NumLevels = numel(levels);

uniformScaleFactor = 1/2 * 1/2; % 2D uniform distribution scale factor

% Pre - Build Flattened Tensor Points
Y_Tensors = cell(NumLevels,1);
W_Tensor = cell(NumLevels,1);

for iLevel = 1:NumLevels
    [yi,wi]=spquad(1,levels(iLevel));

    [Y1, Y2] = meshgrid(yi, yi);

    Y_Tensors{iLevel} = [Y1(:), Y2(:)];

    W_TensorTemp = wi*wi';
    W_Tensor{iLevel} = W_TensorTemp(:);
end

% TODO - not super efficient, could replace with trace(W^T*U(y1,y2))
% OR Reuse function evals
% Evaluate Expected Values

mu_CC = NaN(Nx, NumLevels);
var_CC = NaN(Nx, NumLevels);

KipHalf_vs_x = NaN(Nx,1); % padded with NaNs
KitHalf_vs_x = NaN(Nx,1);

for iLevel = 1:NumLevels   
    NquadPoints = size(Y_Tensors{iLevel},1);
    Iquad_Mu = 0;
    Iquad_uSq = 0;

    Y_CC = Y_Tensors{iLevel};
    W_CC = W_Tensor{iLevel};   

    for i = 1:NquadPoints        
        KitHalf_vs_x(2:end) = getRV_KLE_Fast(dataKLE.x_itHalf, dataKLE, ...
                                    dataKLE.phiMatrix, Y_CC(i,:)'); 
        KipHalf_vs_x(1:end-1) = KitHalf_vs_x(2:end); % Symmetry condition - saves compute
        
        uY = PDE_Solution(uBC, KipHalf_vs_x, KitHalf_vs_x, dx);
   
        Iquad_Mu = Iquad_Mu + uY * W_CC(i);
        Iquad_uSq = Iquad_uSq + uY .^2 * W_CC(i);
    end

    UcrossMoment = Iquad_uSq * uniformScaleFactor;
    mu_CC(:, iLevel) = Iquad_Mu * uniformScaleFactor;
    var_CC(:, iLevel) = UcrossMoment - mu_CC(:, iLevel).^2;
end

plotUxVsNpoints('CC - Stochastic Collocation', x, mu_CC, var_CC, N_levels, false);
plotUxVsNpointsDifference('CC - Stochastic Collocation', x, mu_CC, var_CC, N_levels, false);

%% Test weighting first - e.g. integrate a constant function and see if it scales correctly with my algo





%% Stochastic Collocation - Smolyak


%% Comparison of All Methods



