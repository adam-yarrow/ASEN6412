
% G = gaussian process part of Cmalpha

muG = 0;
targetSigma = 2.1335/10;
sigmaGxx = targetSigma / 5;
lGxx = 40;


tMax = 20;
tMin = 0;

Ktrunc = 3;

t = (tMin:const.dT:tMax)';

aGxx = (tMax - tMin)/2;

meanGFunc = @(t) muG * ones(size(t));

lambdaOptions = struct();
lambdaOptions.mode = 'approximate';
lambdaOptions.k = 200; % Assuming 100 eigenvalues beyond current eval adequate
lambdaOptions.maxIts = 100;
lambdaOptions.fPlot = true;


truncErrorMSE = 0.01; % 1% truncation error

[d, lambda, omega] = findMinNumLambdaKLE(truncErrorMSE, ...
                                    aGxx, lGxx, sigmaGxx, lambdaOptions);


dataKLE = generateKLE_Struct(lambda, omega, meanGFunc, aGxx, d);
[phiMatrix] = MakePhiMatrix(t, dataKLE);

fprintf('min d = %i\n',d);

% Plot realisations
nInstances = 3;
% Y = randn(d, nInstances);

% Truncated Y
Y = trandn(-Ktrunc*ones(d*nInstances,1), Ktrunc*ones(d*nInstances,1));
Y = reshape(Y,d, nInstances);

gd_vs_t = getRV_KLE(t, dataKLE, Y);

plotIndependentRealisations("Part A - l=2", t, gd_vs_t);

%% Test Distributions
nInstances = 1000;
Y = trandn(-Ktrunc*ones(d*nInstances,1), Ktrunc*ones(d*nInstances,1));
Y = reshape(Y,d, nInstances);
gDistrib = getRV_KLE(t, dataKLE, Y);

 figure();
 hold on;
for i = 1:numel(t)

   histogram(gDistrib(i,:));
end


%% Helpers

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
