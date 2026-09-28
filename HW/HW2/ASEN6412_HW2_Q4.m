%{
    TODO - README
%}

rng(8008135); % lol random

%% Common Parameters
% G(x,w) params
sigma = 2;
muG = 1;
meanGFunc = @(t) muG * ones(size(t)); % Constant mean w.r.t. time
varG = sigma^2;

% Support
xMax = 1;
xMin = 0;
dx = 0.001;
x = (xMin:dx:xMax)';
nX = numel(x);

% KLE Options/Parameters
lambdaOptions = struct();
lambdaOptions.mode = 'approximate';
lambdaOptions.k = 10; % Assuming 10 eigenvalues beyond current eval adequate
lambdaOptions.maxIts = 100;
lambdaOptions.fPlot = true;

truncErrorMSE = 0.1; % 10% truncation error allowed
a = 0.5;

%% Part A
% Parameters
l = 2;

% Compute Number of terms for G(x,w) truncation
[d_PartA, lambda_PartA, omega_PartA] = findMinNumLambdaKLE(truncErrorMSE, ...
                                            a, l, sigma, lambdaOptions);
% Plotting Eigenfunctions
nModesToPlot = 4;
fPlotScaledWithLambda = false;
plotEigenFunctions(omega_PartA, lambda_PartA, a, nModesToPlot, fPlotScaledWithLambda)

% Plot 3 independent realisations of G(x,w) from Truncation
nInstances = 3;
Y = randn(d_PartA, nInstances);
dataKLE_PartA = generateKLE_Struct(lambda_PartA, omega_PartA, meanGFunc, a, d_PartA);
gd_vs_x = getRV_KLE(x, dataKLE_PartA, Y);

plotIndependentRealisations("Part A - l=2", x, gd_vs_x);


% Draw samples of G(x,w) and compute sample mean and variance vs x
nInstances = 1E6;
Ylarge = randn(d_PartA, nInstances);
gd_large = getRV_KLE(x, dataKLE_PartA, Ylarge);

meanG_sampled = mean(gd_large,2);
varG_sampled = var(gd_large,0,2);

figure('Name','Sample Moments vs x'); 
subplot(2,1,1);
plot(x, meanG_sampled,'LineWidth',2);
hold on;
yline(muG,'r--','LineWidth',2);
legend('Samples','True Value');
xlabel('x','FontSize',14);
ylabel('E[G]','FontSize',14);
title('Mean of G(x,\omega) vs Position','Interpreter','tex','FontSize',14);
grid on;

subplot(2,1,2);
hold on;
plot(x, varG_sampled,'LineWidth',2);
yline(varG,'r--','LineWidth',2);
legend('Samples','True Value');
xlabel('x','FontSize',14);
ylabel('Var(G)','FontSize',14);
title('Variance of G(x,\omega) vs Position','Interpreter','tex','FontSize',14);
grid on;
ylim([min(varG_sampled), 4.25]);



%% Part B 
% Parameters
l_PartB = 0.2;

% Compute Number of terms for G(x,w) truncation
[d_PartB, lambda_PartB, omega_PartB] = findMinNumLambdaKLE(truncErrorMSE, ...
                                            a, l_PartB, sigma, lambdaOptions);
% Plotting instances
nInstances = 3;
Yb = randn(d_PartB, nInstances);
dataKLE_PartB = generateKLE_Struct(lambda_PartB, omega_PartB, meanGFunc, a, d_PartB);
gd_vs_x_PartB = getRV_KLE(x, dataKLE_PartB, Yb);

plotIndependentRealisations("Part B - l=0.2", x, gd_vs_x_PartB);


%% Helper Functions
function plotIndependentRealisations(partName, x, gd_vs_x)
    figure('Name',sprintf('%s - Plotting Independent Realizations of G',partName));
    hold on;
    for i = 1:size(gd_vs_x,2)    
        plot(x, gd_vs_x(:,i), 'DisplayName',sprintf('\\omega_{%i}', i),'LineWidth',2)
    end
    grid on;
    xlabel('x','FontSize',14)
    ylabel('G(x,\omega)','Interpreter','tex','FontSize',14);
    legend();
    title([partName,' Three Realizations of G(x,\omega) vs x'],'Interpreter','tex','FontSize',14);
end




%% Testing
% Check root finding
omega = 0:0.001:10*2*pi;

omegaSol_odd = [];
omegaSol_Even = [];
omegaSolTanVal_odd = [];
omegaSolTanVal_even = [];
for i = 1:20
    if isEven(i)
        w = getOmegaKLE(i,a,l);
        omegaSol_Even(end+1) = w;
        omegaSolTanVal_even(end+1) = tan(w*a);
    else
        w = getOmegaKLE(i,a,l);
        omegaSol_odd(end+1) = w;
        omegaSolTanVal_odd(end+1) = tan(w*a);
    end
    
end

figure();
subplot(2,1,1);
hold on;
plot(omega/pi, tan(omega*a),'k');
plot(omega/pi, 1./(l*omega), 'r--');
plot(omegaSol_odd/pi, omegaSolTanVal_odd, 'ro');
legend('tan','odd','solutions');
ylim([-2,2]);
grid on;
xlabel('omega (-/\pi)');

subplot(2,1,2);
hold on;
plot(omega/pi, tan(omega*a),'k');
plot(omega/pi, -l*omega,'b--');
plot(omegaSol_Even/pi, omegaSolTanVal_even, 'bo');
legend('tan','even','solutions');
ylim([-100,10]);
grid on;
xlabel('omega (-/\pi)');


%% Testing Lambda i for our problem
maxI = 20;
omegaI = NaN(maxI,1);
lambdaI = NaN(maxI,1);
for i = 1:maxI
    omegaI(i) = getOmegaKLE(i,a,l);
    lambdaI(i) = getLambdaKLE(l,omegaI(i),sigma);
end

figure();
plot(lambdaI);
xlabel('i');
ylabel('lambda_i');
grid on;

%% Testing Minimum d for Lambda
lambdaOptions = struct();
lambdaOptions.mode = 'approximate';
lambdaOptions.k = 10;
% lambdaOptions.MSE = 4;
lambdaOptions.maxIts = 100;
lambdaOptions.fPlot = true;

mseThreshold = 0.1;
[d, lambdaTrunc, omegaTrunc] = findMinNumLambdaKLE(mseThreshold, a, l, sigma, lambdaOptions);


