%{
    This script runs MC simulations for calculating:
        1. grid dependence of mean and variance of u(x) 
        2. Estimates of the mean and variance of u(x)
        3. Probability of having a maximum temperature over specified
        thresholds
        4. Grid dependence of the probability thresholds.

    Note that this code relies on the code that was developed for HW2,
    please see my submission for HW2 for details on the KLE solver etc.
%}
StartupHW3();

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
Nmc = 1E7;
Nx = 101;

% Setup Nmc dependent vars
Y = randn(d, Nmc);
[uBC] = generateU_BC(Nmc, muTheta, sigmaTheta);

x = linspace(xMin,xMax,Nx)';
[u, uStats] = MonteCarloPDE(x, uBC, dataKLE, Y, fPlot);
sigmaEst = sqrt(uStats.var(:,end)); % Get final sigma estimate

% Plotting Convergence
plotReduction('muMidPointReduction', 'varMidPointReduction', {uStats}, Nmc, Nx,...
                'Relative Reduction in <u(0.5)>','Relative Reduction in Var(u(0.5))' );
sgtitle('Relative Reduction in Mean and Variance of u(0.5) vs Number of MC Runs');
plotMeanAndVarForDiffGrid('Mid Point',{uStats},'muAtXhalf','varAtXhalf',...
    '<u(0.5)>','Var(u(0.5))', Nx, Nmc);

% Ploting Final Mean and Variance Estimates
figure('Name','HW3 - u Estimate');
subplot(2,1,1);
hold on;
plot(x, uStats.mu(:,end), 'k','LineWidth',2);
grid on;
xlabel('x');
ylabel('<u(x)> Estimate');
title('MC Estimate of <u(x)>',sprintf('Nmc = %i - Nx = %i', Nmc, Nx))

subplot(2,1,2);
hold on;
plot(x, uStats.var(:,end), 'k','LineWidth',2);
grid on;
xlabel('x');
ylabel('Var(u(x)) Estimate');
title('MC Estimate of Var(u(x))',sprintf('Nmc = %i - Nx = %i', Nmc, Nx))


% Probability of Max Temperature Calculations
% Uses the solution from previous convergence demonstration
[muUmax, varUmax] = CumulativeStats(uStats.max, Nmc);

probCaseA = calculateProbabilityOfMax(muUmax, varUmax, uStats, 3);
probCaseB = calculateProbabilityOfMax(muUmax, varUmax, uStats, 4);


figure('Name','Probability Cases Convergence');
loglog(1:Nmc, probCaseA,...
    'DisplayName','P(u_{max} > <u_{max}> + 3Var(u_{max})^{1/2}');
hold on;
loglog(1:Nmc, probCaseB,...
    'DisplayName','P(u_{max} > <u_{max}> + 4Var(u_{max})^{1/2}');
grid on;
legend();
xlabel('Number of Monte Carlo Runs');
ylabel('Probability of Occurance');
title('Convergence of Probability Estimates Vs Number of MC Runs')

% Reduction Plot Testing
figure();
smoothedReductionProbA = abs(diff(probCaseB)) ./ probCaseB(1:end-1);
loglog(2:Nmc, smoothedReductionProbA)
hold on;
grid on
xlabel('Number of MC Runs')
ylabel('Relative Probability Reduction');
title('Relative Change Between Samples of the Probability')


% Plot Over entire domain
idxToPlot =  [1E1, 1E3, 1E5, 1E7];
figure();

subplot(2,1,1);
hold on;
for i = idxToPlot
    plot(x, uStats.mu(:,i), 'DisplayName',sprintf('N = %i',i));
end
grid on;
xlabel('x');
ylabel('<u(x)> Estimate');
title('MC Estimate of <u(x)>')
legend();

subplot(2,1,2);
hold on;
for i = idxToPlot
    plot(x, uStats.var(:,i), 'DisplayName',sprintf('N = %i',i));
end
grid on;
xlabel('x');
ylabel('Var(u(x)) Estimate');
title('MC Estimate of Var(u(x))')
legend();


%% Grid Convergence Study - On Probabilities
Nmc = 1E5; %% TODO - might need to make this larger
fPlot = false;

% Setup Nmc dependent vars
Y = randn(d, Nmc);
[uBC] = generateU_BC(Nmc, muTheta, sigmaTheta);

% Variables
u = cell(3,1);
uStats = cell(3,1);
x = cell(3,1);
Nx_grid = [5, 21, 81, 101, 161];
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
        if Ngrid == 1
            semilogx(1:Nmc, uStats{i}.(meanVar)(end) + sqrt(uStats{i}.(varVar)(end) ./ (1:Nmc)), ...
                '--k','DisplayName',sprintf('1 sigma Bound'));
            semilogx(1:Nmc, uStats{i}.(meanVar)(end) - sqrt(uStats{i}.(varVar)(end) ./ (1:Nmc)), ...
                '--k','DisplayName',sprintf('1 sigma Bound'));
        end
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

        if Ngrid == 1
            semilogx(2:Nmc, uStats{i}.(varVar)(end) + sqrt(2*uStats{i}.(varVar)(end)^2 ./ (1:Nmc-1)), ...
                '--k','DisplayName',sprintf('1 sigma Bound'));
            semilogx(2:Nmc, uStats{i}.(varVar)(end) -  sqrt(2*uStats{i}.(varVar)(end)^2 ./ (1:Nmc-1)), ...
                '--k','DisplayName',sprintf('1 sigma Bound'));
        end
    end    
    grid on;
    legend();
    xlabel('Number of MonteCarlo Runs');
    ylabel(yVarAxisName,'Interpreter','latex');
    title(['Convergence of ', yVarAxisName, ' vs Nmc Runs'],'Interpreter','latex')

    sgtitle(sprintf('Convergence of %s vs Number of MC Runs',figName));
end




function plotReduction(varNameMu, varNameVariance, uStats, Nmc, Nx, yLabelMu, yLabelVar)
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
    ylabel(yLabelMu);

    subplot(2,1,2);
    for i = 1:Ngrid
        semilogy(1:Nmc, uStats{i}.convergence.(varNameVariance),'DisplayName',sprintf('Nx = %i',Nx(i)));
        hold on;
    end
    grid on;
    legend();
    xlabel('Number of MonteCarlo Runs')
    ylabel(yLabelVar);
end
