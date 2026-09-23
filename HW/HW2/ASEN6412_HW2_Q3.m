x1 = 1;
x2 = 0:0.001:10;

Fx2_Given_X1 = 1 - (1+x2).*exp(-x2.*(1+x1));

figure();
plot(x2, Fx2_Given_X1);


%% Testing domain of Wk
x1 = 0:0.1:1;
dz2 = 0.001;
z2 = 0:dz2:1-dz2;

[X1, Z2] = meshgrid(x1,z2);

a = log(1-Z2) - (1+X1);
Z = -(1+X1).*exp(a);

figure();
plot3(X1,Z2,Z);
xlabel('X1')
ylabel('Z2');
grid on;