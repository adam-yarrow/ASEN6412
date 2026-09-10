%% Question 4 - Part A testing
dx = 0.01;
dy = 0.01;

x = -2:dx:2;
nX = numel(x);
y = -1.5:dy:1.5;
nY = numel(y);


[X, Y] = meshgrid(x,y);

fXY = zeros(size(X));

validIdxX = 1:nX;
validIdxX = validIdxX((x >= 0) & (x <= 1));
validIdxY = 1:nY;
validIdxY = validIdxY((y >= 0) & (y <= 1));

for i = 1:numel(validIdxY)
    for j = 1:numel(validIdxX)
        fXY(validIdxY(i), validIdxX(j)) = 2*(1-x(validIdxX(j)));
    end
end

% pdf
figure();
% plot3(X,Y,fXY,'r*');
surf(X,Y,fXY);
xlabel('x')
ylabel('y')
zlabel('fXY');
grid on;
title('pdf');



%% Part d
dz = 0.01;
Z = -1:dz:3;
nZ = numel(Z);

fZ = zeros(nZ,1);

for i = 1:nZ
    z = Z(i);
    if z <= 1 & z >= 0
        fZ(i) = z*(2-z);
    elseif z > 1 & z <=2
        fZ(i) = (z-2)^2;
    end    
end

FZ = cumsum(fZ*dz);


% Checking pdf
Npoints = 10^5;
Y_RV = rand(Npoints,1); % Uniform random
xRV_draws = rand(Npoints,1);
X_RV = 1 - (1 - xRV_draws).^0.5; % inverse CDF distribution

Z_True = Y_RV + X_RV;

figure('Name','Q4 - Part D');
subplot(1,2,1);
plot(Z, fZ);
hold on;
histogram(Z_True,'Normalization','pdf');
title('pdf');

subplot(1,2,2);
plot(Z, FZ);
title('cdf');



