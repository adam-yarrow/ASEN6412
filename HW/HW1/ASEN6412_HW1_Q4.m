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


% CDF
