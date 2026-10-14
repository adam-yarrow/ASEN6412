function [u, uStats] = MonteCarloPDE(x, uBC, dataKLE, Y, fPlot)
%{
    Function to evaluation the MC approximation of the PDE for a fixed grid
    size.

    x = grid to evaluate over
    uBC = Nmc x 2 vector containing the boundary conditions for each run
    dataKLE = KLE data structure
    Y = random normal variables used for KLE expansion of G(x,w), size d x
    Nmc
    fPlot = Flag to enable/disable plotting
%}
dx = x(2) - x(1);
Nx = numel(x);
Nmc = size(uBC,1);

%% Sampling K(x)
x_itHalf = x - dx/2; 
x_itHalf = x_itHalf(2:end); % First entry outside of domain -> remove
KitHalf_vs_x = exp(getRV_KLE(x_itHalf, dataKLE, Y)); % Exponential K(x,w) = exp(Gd(x,w))
KipHalf_vs_x = KitHalf_vs_x; % Symmetry condition - saves compute

% pad with Nans to make size Nx1
KipHalf_vs_x = [KipHalf_vs_x; NaN(1,Nmc)]; 
KitHalf_vs_x = [NaN(1,Nmc); KitHalf_vs_x];

%% Solve MC
u = NaN(Nx,Nmc);
for i = 1:Nmc
    u(:,i) = HarmonicPDE_Solution(uBC(i,:), ...
                KipHalf_vs_x(:,i), KitHalf_vs_x(:,i), dx);
end

%% Generate Statistics
[muU, varU] = CumulativeStats(u, Nmc);

% max u(x)
umax = max(u,[],1);

%% Convergence Criteria
% Convergence on mean and variance of u(x) at x = 0.5;
xHalfIdx = find(x == 0.5); %% Brittle - fix
muU_xHalfReduction = calculateReduction(muU(xHalfIdx, :));
varU_xHalfReduction = calculateReduction(varU(xHalfIdx, :));

% Convergence on integral of mu(x) and var(x)
integralMuU = trapz(x,muU,1);
integralVarU = trapz(x,varU,1);

% Integral of difference between different number of runs


muU_IntegralReduction = calculateReduction(integralMuU);
varU_IntegralReduction = calculateReduction(integralVarU);



%% Logging Data
uStats.mu = muU;
uStats.var = varU;
uStats.max = umax;
uStats.convergence.muMidPointReduction = muU_xHalfReduction;
uStats.convergence.varMidPointReduction = varU_xHalfReduction;
uStats.convergence.muIntegralReduction = muU_IntegralReduction;
uStats.convergence.varIntegralReduction = varU_IntegralReduction;


%% Plotting Convergence
if fPlot


    %% TODO - add plots for convergence criteria here



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
end

end

function reductionValue = calculateReduction(var)    
    % reductionValue = abs([NaN, diff(var)]' ./ var');
    reductionValue = abs([NaN, diff(var)]');
end