

% G = gaussian process part of Cmalpha


tMax = 10;
tMin = 0;


muG = 0;
targetSigma = 2.1335/10;
sigmaGxx = targetSigma * 0.75; %/5
lGxx = tMax*2*0.32;


Ktrunc = 1; %3;

t = (tMin:const.dT:tMax)';

aGxx = (tMax - tMin)/2;

meanGFunc = @(t) muG * ones(size(t));

lambdaOptions = struct();
lambdaOptions.mode = 'approximate';
lambdaOptions.k = 200; % Assuming 100 eigenvalues beyond current eval adequate
lambdaOptions.maxIts = 100;
lambdaOptions.fPlot = false;


truncErrorMSE = 0.1; % 10% truncation error

[d, lambda, omega] = findMinNumLambdaKLE(truncErrorMSE, ...
                                    aGxx, lGxx, sigmaGxx, lambdaOptions);

lambda = lambda(1:d);
omega = omega(1:d);
dataKLE = generateKLE_Struct(lambda, omega, meanGFunc, aGxx, d);
[phiMatrix] = MakePhiMatrix(t, dataKLE);

fprintf('min d = %i\n',d);

% Plot realisations
nInstances = 10;
% Y = randn(d, nInstances);

% Truncated Y
% Y = trandn(-Ktrunc*ones(d*nInstances,1), Ktrunc*ones(d*nInstances,1));
% Y = reshape(Y,d,nInstances);

Y = rand(d,nInstances)*2 - 1; % trying uniform Yi = [-1,1]

gd_vs_t = NaN(numel(t), nInstances);
for i = 1:nInstances
    gd_vs_t(:,i) = getRV_KLE(t, dataKLE, Y(:,i));
end

plotIndependentRealisations("Part A - l=2", t, gd_vs_t);


figure();
plot(t', (phiMatrix.*lambda(1:d)));
grid on;

%% Test Distributions
nInstances = 1000;
Y = trandn(-Ktrunc*ones(d*nInstances,1), Ktrunc*ones(d*nInstances,1));
Y = reshape(Y,d, nInstances);
gDistrib = getRV_KLE(t, dataKLE, Y);

 figure();
 hold on;
 % cnt = 0
for i = 1:numel(t)
    if mod(i,50)
        histogram(gDistrib(i,:));
        % cnt = cnt + 1
    end
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
