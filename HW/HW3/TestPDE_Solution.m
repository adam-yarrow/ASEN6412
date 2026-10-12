%{
Tests the PDE solver for constant K.
%}
%% General Configuration
uBC = [0,0];

dx = 0.01;
x = (0:dx:1)';
N = numel(x);

Kconst = 5;

%% Case 1: K(x) = Constant
K = ones(N,1)*Kconst;

uAnalytical_ConstK = ConstantKSolution(Kconst,x);
uDiscrete_ConstK = HarmonicPDE_Solution(uBC, K, dx);
percErr_ConstK = (uDiscrete_ConstK - uAnalytical_ConstK) ./ uAnalytical_ConstK * 100;

plotComparison(x, uDiscrete_ConstK, uAnalytical_ConstK, percErr_ConstK, 'Constant K')

%% Case 2: K(x) = K*(1+x)
K = Kconst*(ones(N,1) + x);

uAnalytical_LinK = LinearKSolution(Kconst,x);
uDiscrete_LinK= HarmonicPDE_Solution(uBC, K, dx);
percErr_LinK = (uDiscrete_LinK - uAnalytical_LinK) ./ uAnalytical_LinK * 100;

plotComparison(x, uDiscrete_LinK, uAnalytical_LinK, percErr_LinK, 'Linear K')



%% Analytical Solutions
function u = ConstantKSolution(K,x)
    % Assuming K(x) = K;
    u = (x/(2*K)) .* (1-x);    
end

function u = LinearKSolution(K,x)
    % Assuming K(x) = K*(1+x)
    cp1 = 1/log(2);
    u = (1/K)*(cp1*log(1+x) - x);
end

%% Helper Functions
function plotComparison(x, uDiscrete, uAnalytical, percErr, Ktype)
    figure('Name',sprintf('HW3 - PDE Verification - %s',Ktype)); 
    subplot(2,1,1);
    hold on;
    plot(x, uDiscrete,'b*');
    plot(x, uAnalytical,'r');
    grid on;
    legend('Discrete','Analytical');
    title(sprintf('PDE %s Solution Verification - Values',Ktype));
    xlabel('x');
    ylabel('u(x)');
    
    subplot(2,1,2);
    plot(x, percErr);
    grid on;
    title(sprintf('PDE %s Solution Verification - Error',Ktype));
    xlabel('x');
    ylabel('Error (%)');
end