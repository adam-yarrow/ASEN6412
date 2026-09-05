%% Testing Part A
X = randn(1E8,1);

Y = (1/sqrt(2)) * (X.^2 - 1);

yMin = -1/sqrt(2);
yMean = mean(Y);
yVar = var(Y);


% Test distribution function
y = -0.70:0.01:25;

pdfY = exp(1).^(-0.5*(sqrt(2)*y+1)) ./ (sqrt(pi)*(sqrt(2)*y+1).^0.5);

figure(); 
histogram(Y, 'Normalization','pdf');
hold on;
plot(y, pdfY);

%% Testing Part D
rhoXY = corrcoef(X,Y);


