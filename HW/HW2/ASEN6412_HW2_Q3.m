%{
    This script generates samples of X1 and X2 according to the joint
    probability distribution function for Q3 of HW2. It then plots the
    histograms of the samples and compares to the true CDF and pdf results.
%}
rng(121246554);

%% Generating Samples Via Inverse Method
Nsamples = 1E6; % Number of samples to draw

% Generate IID Uniform z = [z1, z2]^T
z = rand(2,Nsamples);

% Generate X1 from inverse CDF
x1 = Finv_X1(z(1,:));

% Generate X2 from conditional inverse CDF
x2 = Finv_X2givenX1(x1, z(2,:));

% Truth Joint PDF
dx = 0.1;
x1_true = 0:dx:10;
x2_true = 0:dx:10;
[X1_true, X2_true] = meshgrid(x1_true, x2_true);

% Analytical Solution for E[X1X2]
E_X1X2_Analytical = exp(1)*expint(1);

% Estimate density of samples
x = [x1', x2'];
xi = [X1_true(:), X2_true(:)];
[fXEst_Hist, xEdges, yEdges] = histcounts2(x1,x2,'Normalization','pdf');
xCentres = xEdges(1:end-1) + (xEdges(2:end) - xEdges(1:end-1))./2;
yCentres = yEdges(1:end-1) + (yEdges(2:end) - yEdges(1:end-1))./2;
[XCentres, YCentres] = meshgrid(xCentres,yCentres);

% Statistics Extraction
mu1 = mean(x1);
mu2 = mean(x2);

cov12 = cov(x1,x2);
E_X1X2 = cov12(1,2) + mu1*mu2; % Adding product of means back to get cross moment


%% Plotting
% Plotting - Density
figure('Name','HW2 - Q3 - Joint PDF');
hold on;
surf(X1_true, X2_true, f_X2_X1(X1_true, X2_true),'EdgeColor','none');
plot3(XCentres, YCentres, fXEst_Hist,'r.');
grid on;
legend('Analytical','Inverse Method - Histogram')
xlabel('X1');
ylabel('X2');
zlabel('f_{X1,X2}');
title('Joint pdf - Analytical vs Sampled');
xlim([min(x1_true), max(x1_true)]);
ylim([min(x2_true), max(x2_true)]);

% Plotting Moments vs Time
n = 1:Nsamples;
mu1_vsN = cumsum(x1,2) ./ n;
mu2_vsN = cumsum(x2,2) ./ n;
E_X1X2_vsN = cumsum(x1.*x2,2) ./ n;

figure('Name','HW2 - Q3 - Moments vs Number of Samples')
subplot(2,1,1);
semilogx(n, mu1_vsN,'r');
hold on;
semilogx(n, mu2_vsN,'b');
semilogx([1, Nsamples],[1 1],'k--');
legend('\mu_1','\mu_2','Analytical Solution');
xlabel('Number of Samples');
ylabel('Running Mean');
grid on;
title('Mean vs Number of Samples for X_1 and X_2')

subplot(2,1,2);
semilogx(n, E_X1X2_vsN,'r');
hold on;
semilogx([1, Nsamples],[E_X1X2_Analytical E_X1X2_Analytical],'k--');
legend('E[X_1X_2]','Analytical Solution');
xlabel('Number of Samples');
ylabel('Running E[X_1X_2]');
grid on;
title('E[X_1X_2] vs Number of Samples')



%% Helper Functions
% Marginal CDF of X1
function c = F_X1(x1)
    c = 1 - exp(-x1);
end

% Marginal pdf of X1
function p = f_X1(x1)
    p = exp(-x1);
end

% Conditional CDF of X2 given X1
function c = F_X2givenX1(x1,x2)
    c = 1 - (1+x2).*exp(-x2.*(1+x1));
end

% Joint PDF
function p = f_X2_X1(x1,x2)
    p = (x1 + x2 + x1.*x2) .* exp(-x1 - x2 - x1.*x2);
end

% Inverse CDF of X1
function x1 = Finv_X1(z1)
    x1 = -log(1-z1);
end

% Inverse CDF of X2 given X1
function x2 = Finv_X2givenX1(x1,z2)
    a = log(1-z2) - (1+x1);

    % Lambert Soln
    k = -1; % Branch
    z = -(1+x1).*exp(a);
    W = lambertw(k, z);

    x2 = exp(a - W) - 1;    
end
