x1 = 1;
x2 = 0:0.001:10;

Fx2_Given_X1 = 1 - (1+x2).*exp(-x2.*(1+x1));

figure();
plot(x2, Fx2_Given_X1);