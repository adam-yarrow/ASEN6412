%{
    TODO - update README
%}
% Fix random seed for repeatability
rng(12102026);

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

%% Grid Convergence Study - Mean and Variance
Nmc = 1E4;

fPlot = false;

% Setup Nmc dependent vars
Y = randn(d, Nmc);
[uBC] = generateU_BC(Nmc, muTheta, sigmaTheta);

% Variables
maxGrids = 5;
u = cell(maxGrids,1);
uStats = cell(maxGrids,1);
x = cell(maxGrids,1);
Nx = NaN(maxGrids,1);
for i = 1:maxGrids
    Nx(i) = 4^i + 1;
    x{i} = linspace(xMin,xMax,Nx(i))';
    [u{i}, uStats{i}] = MonteCarloPDE(x{i}, uBC, dataKLE, Y, fPlot);
end

% plotReduction('muMidPointReduction', 'varMidPointReduction', uStats, Nmc, Nx);
% plotReduction('muIntegralReduction', 'varIntegralReduction', uStats, Nmc, Nx);
plotMeanAndVarForDiffGrid('Mid Point',uStats,'muAtXhalf','varAtXhalf',...
    '<u(0.5)>','Var(u(0.5))', Nx, Nmc);

% plotMeanAndVarForDiffGrid('Spatial Mean',uStats,'integralMu','integralVarU',...
%     '$\int_{0}^{1}<u(x)>dx$','$\int_{0}^{1}Var(u(x))dx$', Nx, Nmc);

%% Convergence of Mean and Variance of u(x)
Nmc = 1E6;
Nx = 101;

% Setup Nmc dependent vars
Y = randn(d, Nmc);
[uBC] = generateU_BC(Nmc, muTheta, sigmaTheta);

x = linspace(xMin,xMax,Nx)';
profile on
[u, uStats] = MonteCarloPDE(x, uBC, dataKLE, Y, fPlot);
profile viewer
sigmaEst = sqrt(uStats.var(:,end));

figure('Name','Convergence of Mean and Variance');
%TODO - put convergence criteria here...

% Ploting u(x)
figure('Name','HW3 - u Estimate');
hold on;
plot(x, uStats.mu(:,end), 'r');
plot(x, uStats.mu(:,end) + sigmaEst,'r--');
plot(x, uStats.mu(:,end) - sigmaEst,'r--');
legend('\mu','\mu +/- \sigma');
grid on;
xlabel('x');
ylabel('u(x) Estimate');
title(sprintf('Monte Carlo Estimate of u(x) - Mean and Variance - N = %i', Nmc))


% Probability of Max Temperature Calculations
% Uses the solution from previous convergence demonstration
[muUmax, varUmax] = CumulativeStats(uStats.max, Nmc);

probCaseA = calculateProbabilityOfMax(muUmax, varUmax, uStats, 3);
probCaseB = calculateProbabilityOfMax(muUmax, varUmax, uStats, 4);


figure('Name','Probability Cases Convergence');
semilogx(1:Nmc, probCaseA,...
    'DisplayName','P(u_{max} > <u_{max}> + 3Var(u_{max})^{1/2}');
hold on;
semilogx(1:Nmc, probCaseB,...
    'DisplayName','P(u_{max} > <u_{max}> + 4Var(u_{max})^{1/2}');
grid on;
legend();
xlabel('Number of Monte Carlo Runs');
ylabel('Probability of Occurance');
title('Convergence of Probability Estimates Vs Number of MC Runs')


profile off


%% Testing prob convergence
figure();
smoothedReduction = movmean(abs(diff(probCaseA)) ./ probCaseA(1:end-1),10000);
semilogy(2:Nmc, smoothedReduction)
grid on
xlabel('Number of MC Runs')
ylabel('Smoothed Difference');
title('Smoothed Absolute Change Between Samples of the Probability')


%% Grid Convergence Study - On Probabilities
Nmc = 1E6; %% TODO - might need to make this larger
fPlot = false;

% Setup Nmc dependent vars
Y = randn(d, Nmc);
[uBC] = generateU_BC(Nmc, muTheta, sigmaTheta);

% Variables
u = cell(3,1);
uStats = cell(3,1);
x = cell(3,1);
Nx_grid = [5, 21, 81, 161];%[11, 51, 101, 501, 1001];
Ngrid = numel(Nx_grid);
probCaseA_Grid = NaN(Nmc, Ngrid);
probCaseB_Grid = NaN(Nmc, Ngrid);


for i = 1:Ngrid
    x{i} = linspace(xMin,xMax,Nx_grid(i))';
    [~, uStats{i}] = MonteCarloPDE(x{i}, uBC, dataKLE, Y, fPlot);

    [muUmax, varUmax] = CumulativeStats(uStats{i}.max, Nmc);

    probCaseA_Grid(:,i) = calculateProbabilityOfMax(muUmax, varUmax, uStats{i}, 3);
    probCaseB_Grid(:,i) = calculateProbabilityOfMax(muUmax, varUmax, uStats{i}, 4);
end

figure('Name','Grid Refinement Study - Probabilities'); 
subplot(2,1,1);
for i = 1:Ngrid
    loglog(1:Nmc, probCaseA_Grid(:,i),'DisplayName',sprintf('Nx = %i',Nx_grid(i)));
    hold on;
end
grid on;
legend();
xlabel('Number of MonteCarlo Runs')
ylabel('Pr(u_{max} > <u_{max}> + 3\sigma_{u_{max}}');

subplot(2,1,2);
for i = 1:Ngrid
    loglog(1:Nmc, probCaseB_Grid(:,i),'DisplayName',sprintf('Nx = %i',Nx_grid(i)));
    hold on;
end
grid on;
legend();
xlabel('Number of MonteCarlo Runs')
ylabel('Pr(u_{max} > <u_{max}> + 4\sigma_{u_{max}}');
sgtitle('Grid Convergence of Probability of Maximum');


% Reduction Plot
figure();
for i = 1:Ngrid
    smoothedReductionProbA = movmean(abs(diff(probCaseA_Grid(:,i))) ./ probCaseA_Grid(1:end-1,i),10000);
    semilogy(2:Nmc, smoothedReductionProbA)
    hold on;
end
grid on
xlabel('Number of MC Runs')
ylabel('Smoothed Difference');
title('Smoothed Absolute Change Between Samples of the Probability')



%% Plotting Functions

function plotMeanAndVarForDiffGrid(figName,uStats,meanVar,varVar,yMuAxisName,yVarAxisName,Nx,Nmc)
    Ngrid = numel(Nx);
    figure('Name',sprintf('Convergence of %s vs Number of MC Runs',figName)); 
    subplot(2,1,1);
    for i = 1:Ngrid
        semilogx(1:Nmc, uStats{i}.(meanVar),'DisplayName',sprintf('Nx = %i',Nx(i)));
        hold on;
    end
    grid on;
    legend();
    xlabel('Number of MonteCarlo Runs')
    ylabel(yMuAxisName,'Interpreter','latex');
    title(['Convergence of ', yMuAxisName, ' vs Nmc Runs'],'Interpreter','latex')

    subplot(2,1,2);
    for i = 1:Ngrid
        semilogx(1:Nmc, uStats{i}.(varVar),'DisplayName',sprintf('Nx = %i',Nx(i)));
        hold on;
    end
    grid on;
    legend();
    xlabel('Number of MonteCarlo Runs');
    ylabel(yVarAxisName,'Interpreter','latex');
    title(['Convergence of ', yVarAxisName, ' vs Nmc Runs'],'Interpreter','latex')

    sgtitle(sprintf('Convergence of %s vs Number of MC Runs',figName));
end




function plotReduction(varNameMu, varNameVariance, uStats, Nmc, Nx)
    Ngrid = numel(Nx);
    figure(); 
    subplot(2,1,1);
    for i = 1:Ngrid
        semilogy(1:Nmc, uStats{i}.convergence.(varNameMu),'DisplayName',sprintf('Nx = %i',Nx(i)));
        hold on;
    end
    grid on;
    legend();
    xlabel('Number of MonteCarlo Runs')
    % ylabel('$|(\bar{u_{i}}-\bar{u_{i-1}})/\bar{u_{i-1}}|$','Interpreter','latex')
    subplot(2,1,2);
    for i = 1:Ngrid
        semilogy(1:Nmc, uStats{i}.convergence.(varNameVariance),'DisplayName',sprintf('Nx = %i',Nx(i)));
        hold on;
    end
    grid on;
    legend();
    xlabel('Number of MonteCarlo Runs')
end
