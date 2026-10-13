%{
    This script generates monte-carlo samples of the stochastic PDE
    d/dx(K(x,w)*du(x,w)/dx) = 1.0
 
    K(x,w) = exp(G(x,w)
    x = [0,1]
    G(x,w) = Gaussian random process with mean of 1.0 and exponential
    covariance function.

    TODO - might make this a function later? If its worth it/can be arsed
%}
% Fix random seed for repeatability
rng(12102026);

%% Variables
Nx = 11; % Number of x points
Nmc = 1; 1E3; %% TODO bump these numbers up

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

% x = 0 BC
uBC_0 = 0;

x =  linspace(0,1,Nx)';
dx = x(2) - x(1);

%% Setup Stochastic BC
uBC_1 = muTheta + sigmaTheta*randn(Nmc,1);

%% Generate K(x,w) samples
% Get d eigenfunctions
lambdaGxx = NaN(d,1);
omegaGxx = NaN(d,1);
aGxx = (x(end) - x(1))/2; % Calculate domain boundary

for i = 1:d
    omegaGxx(i) = getOmegaKLE(i,aGxx,lGxx);
    lambdaGxx(i) = getLambdaKLE(lGxx,omegaGxx(i),sigmaGxx);
end

% Generate Random Samples
Y = randn(d, Nmc);
meanGFunc = @(x) muG * ones(size(x));
dataKLE = generateKLE_Struct(lambdaGxx, omegaGxx, meanGFunc, aGxx, d);

% x_ipHalf = x + dx/2; x_ipHalf = x_ipHalf(1:end-1);
x_itHalf = x - dx/2; x_itHalf = x_itHalf(2:end);
% KipHalf_vs_x = exp(getRV_KLE(x_ipHalf, dataKLE, Y));
KitHalf_vs_x = exp(getRV_KLE(x_itHalf, dataKLE, Y));
KipHalf_vs_x = KitHalf_vs_x; % Symmetry condition - saves compute

% pad with Nans to make size Nx1
KipHalf_vs_x = [KipHalf_vs_x; NaN(1,Nmc)]; 
KitHalf_vs_x = [NaN(1,Nmc); KitHalf_vs_x];

%% Solve MC
u = NaN(Nx,Nmc);
for i = 1:Nmc
    uBC = [uBC_0, uBC_1(i)];
    u(:,i) = HarmonicPDE_Solution(uBC,KipHalf_vs_x(:,i), KitHalf_vs_x(:,i),dx);
end

%% Generate Statistics
% Running Mean
muU = cumsum(u,2) ./ (1:Nmc);
% Running Var
uSquared = cumsum(u.^2,2);
NmeanU = (1:Nmc) .* muU.^2;
varU = (uSquared - NmeanU) ./ ((1:Nmc) - 1);
varU(:,1) = NaN;

%% TODO - calculate SEM of muU, maybe any stats on varU?



%% Plotting Convergence
[X, N_MC] = meshgrid(x, 1:Nmc);

figure();
% surf(X,N_MC,muU','EdgeAlpha',0.1);
plot3(X,N_MC,muU');
grid on;
xlabel('x');
ylabel('Number of Monte Carlo Samples');
zlabel('E[u(x)]');
title('Mean of u(x) vs x and Number of Samples')


figure();
% surf(X,N_MC,varU','EdgeAlpha',0.1);
plot3(X,N_MC,varU');
grid on;
xlabel('x');
ylabel('Number of Monte Carlo Samples');
zlabel('Var(u(x))');
title('Variance of u(x) vs x and Number of Samples')

% plot 1-sigma on meanU and varU estimates?

%% TODO - maybe pick a few points along x and plot these individually with SEM bounds?

%% Probability of Max Temperature Calculations