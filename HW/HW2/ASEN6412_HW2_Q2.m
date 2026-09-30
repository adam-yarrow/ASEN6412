%{
    HW2 - Q2
    This script generates samples of the Chebyshev probability density
    using the inversion method.
%}

%% Inputs
nSamples = 10000;

%% Inversion Method
%% 
rng(100); % Fix seed for repeatability
u = rand(nSamples,1);
xInv = Finv(u);

%% Estimate Density and distributions functions
xGridEst = -1.0:0.001:1.0;
[fXEst, xEst] = ksdensity(xInv,xGridEst,'Bandwidth',0.001); %'Support',[-1,1],
%% TODO - play with bandwidth for smoothing and potentially the support?

[FXEst, xEst_FX] = ecdf(xInv);

% Generate Real CDF and PDF
xGrid = -1.5:0.001:1.5;
pdfTrue = fX(xGrid);
cdfTrue = FX(xGrid);


%% Plotting
figure('Name','HW2 - Q2 - PDF'); 
hold on;
plot(xGrid, pdfTrue,'LineWidth',2);
plot(xEst, fXEst,'LineWidth',1);
grid on;
xlabel('x','FontSize',14);
ylabel('f_X(x)','FontSize',14);
legend('True pdf','Inverse Method','FontSize',14)
title('pdf - Drawing from pdf via Inversion Method','FontSize',14)

figure('Name','HW2 - Q2 - CDF'); 
hold on;
plot(xGrid, cdfTrue,'LineWidth',2);
plot(xEst_FX, FXEst,'LineWidth',2);
grid on;
xlabel('x','FontSize',14);
ylabel('F_X(x)','FontSize',14);
legend('True CDF','Inverse Method','FontSize',14)
title('CDF - Drawing from pdf via Inversion Method','FontSize',14)


%% Helpers
% Inverse CDF
function x = Finv(u)
    % Protection never hurt anyone
    if u > 1
        u = 1;
    elseif u < 0
        u = 0;
    end

    x = sin(pi.*(u + 0.5));    
end

% True PDF
function p = fX(x)
    p = zeros(size(x));

    % Density in -1 < x < 1
    validXidx = (x < 1) & (x > -1);
    validX = x(validXidx);

    p(validXidx) = 1 ./ (pi * (1 - validX.^2).^0.5); 
end

% True CDF

function c = FX(x)
    c = zeros(size(x));

    % CDF in -1<= x <= 1
    validXidx = (x < 1) & (x > -1);
    validX = x(validXidx);

    c(validXidx) = asin(validX)/pi + 0.5;

    % CDF in x >= 1
    c(x >= 1) = 1;
end