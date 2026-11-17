% Austin Bodin
% Adam Yarrow

addpath(genpath(pwd));
addpath(genpath('..\..\Library'));


load("u_samples.mat");
load("y_samples.mat");

p = 3;
d = 20;
P = factorial(p+d) / (factorial(p) * factorial(d)) - 1;
Ncases = ceil([0.2; 0.5; 0.9] * (P+1));
numN = numel(Ncases);

opts = spgSetParms('iterations',6000,'verbosity',1,...
     'optTol',1e-9,'bpTol',1e-9);

sigma = 0.005;

index_pc = nD_polynomial_array(d,p);
c_ell1 = cell(numN,1);

meanU = NaN(numN,1);
varU = NaN(numN,1);

for iN = 1:numN
    N = Ncases(iN);

    % Make psi
    psi = zeros(N, P+1);
    for i = 1:N
        psi(i,:) = piset(y_samples(i,:), index_pc);
    end

    delta = sigma*norm(u_samples(1:N));

    % solve system   
    [c_ell1{iN},res,gap,info] = spg_bpdn(psi,u_samples(1:N),delta,opts);

    % extract Mean and variance
    meanU(iN) = c_ell1{iN}(1); % First order coefficient

    varU(iN) = sum(c_ell1{iN}(2:end).^2); % Remove first term in summation to remove mean^2
    
end

%% Plotting

figure();
hold on;
plot(Ncases, meanU);
grid on;
xlabel('Number of Samples');
ylabel('E[u(y)]');
title('Mean of u(y) vs Number of Samples - l_1 Minimization')


figure();
hold on;
plot(Ncases, varU);
grid on;
xlabel('Number of Samples');
ylabel('Var(u(y))');
title('Variance of u(y) vs Number of Samples - l_1 Minimization')