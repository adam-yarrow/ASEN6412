function plotEigenFunctions(omega, lambda, a, nToPlot, fPlotScaled)
xTilde = linspace(-a,a,1000);

if fPlotScaled
    scaling = sqrt(lambda);
    scalingType = 'Lambda Scaled ';
else
    scaling = ones(size(lambda));
    scalingType = '';
end

figure('Name','Plotting Eigenfunctions');
hold on;
for i = 1:nToPlot    
    phiFunc = evalEigenfunctionKLE(omega(i), a, i);
    plot(xTilde, scaling(i)*phiFunc(xTilde),...
            'DisplayName',sprintf('i=%i',i));
end
legend();
xlabel('x (-a<x<a)');
ylabel('\phi_i(x)');
grid on;
title(sprintf('%sEigenfunctions vs non-dimensional x',scalingType));
end