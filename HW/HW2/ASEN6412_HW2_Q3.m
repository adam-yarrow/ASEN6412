% x1 = 1;
% x2 = 0:0.001:10;
% 
% Fx2_Given_X1 = 1 - (1+x2).*exp(-x2.*(1+x1));
% 
% figure();
% plot(x2, Fx2_Given_X1);
% 
% 
% %% Testing domain of Wk
% x1 = 0:0.1:1;
% dz2 = 0.001;
% z2 = 0:dz2:1-dz2;
% 
% [X1, Z2] = meshgrid(x1,z2);
% 
% a = log(1-Z2) - (1+X1);
% Z = -(1+X1).*exp(a);
% 
% figure();
% plot3(X1,Z2,Z);
% xlabel('X1')
% ylabel('Z2');
% grid on;


%{
    This script generates samples of X1 and X2 according to the joint
    probability distribution function for Q3 of HW2. It then plots the
    histograms of the samples and compares to the true CDF and pdf results.
%}
rng(121246554)

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


% Estimate density of samples
x = [x1', x2'];
xi = [X1_true(:), X2_true(:)];
[fXEst_Hist, xEdges, yEdges] = histcounts2(x1,x2,'Normalization','pdf');
xCentres = xEdges(1:end-1) + (xEdges(2:end) - xEdges(1:end-1))./2;
yCentres = yEdges(1:end-1) + (yEdges(2:end) - yEdges(1:end-1))./2;
[XCentres, YCentres] = meshgrid(xCentres,yCentres);
% fXEst = ksdensity(x,xi,'Bandwidth',0.1); %% Kinda works but is bad

% Plotting
figure('Name','HW2 - Q3 - Joint PDF');
hold on;
surf(X1_true, X2_true, f_X2_X1(X1_true, X2_true),'EdgeColor','none');
% plot3(xi(:,1), xi(:,2), fXEst,'r*');
plot3(XCentres, YCentres, fXEst_Hist,'r*');
% histogram2(x1,x2,'Normalization','pdf')
grid on;
legend('Analytical','Inverse Method - Histogram')
xlabel('X1');
ylabel('X2');
zlabel('f_{X1,X2}');
title('Joint pdf - Analytical vs Sampled');
xlim([min(x1_true), max(x1_true)]);
ylim([min(x2_true), max(x2_true)]);


% Statistic Extraction
mu1 = mean(x1);
mu2 = mean(x2);

cov12 = cov(x1,x2);
E_X1X2 = cov12(1,2) + mu1*mu2; % Adding product of means back to get cross moment


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
