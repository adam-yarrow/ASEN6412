%{
    TODO - README
%}

%% Parameters - Part A
a = 0.5;
l = 2;
sigma = 2;
truncError = 0.1;

%% Testing
% Check root finding
omega = 0:0.001:10*2*pi;

omegaSol_odd = [];
omegaSol_Even = [];
omegaSolTanVal_odd = [];
omegaSolTanVal_even = [];
for i = 1:20
    if isEven(i)
        w = getOmegaKLE(i,a,l);
        omegaSol_Even(end+1) = w;
        omegaSolTanVal_even(end+1) = tan(w*a);
    else
        w = getOmegaKLE(i,a,l);
        omegaSol_odd(end+1) = w;
        omegaSolTanVal_odd(end+1) = tan(w*a);
    end
    
end

figure();
subplot(2,1,1);
hold on;
plot(omega/pi, tan(omega*a),'k');
plot(omega/pi, 1./(l*omega), 'r--');
plot(omegaSol_odd/pi, omegaSolTanVal_odd, 'ro');
legend('tan','odd','solutions');
ylim([-2,2]);
grid on;
xlabel('omega (-/\pi)');

subplot(2,1,2);
hold on;
plot(omega/pi, tan(omega*a),'k');
plot(omega/pi, -l*omega,'b--');
plot(omegaSol_Even/pi, omegaSolTanVal_even, 'bo');
legend('tan','even','solutions');
ylim([-100,10]);
grid on;
xlabel('omega (-/\pi)');


%% Testing Lambda i for our problem
maxI = 20;
omegaI = NaN(maxI,1);
lambdaI = NaN(maxI,1);
for i = 1:maxI
    omegaI(i) = getOmegaKLE(i,a,l);
    lambdaI(i) = getLambdaKLE(l,omegaI(i),sigma);
end

figure();
plot(lambdaI);
xlabel('i');
ylabel('lambda_i');
grid on;

%% Testing Minimum d for Lambda
lambdaOptions = struct();
lambdaOptions.mode = 'approximate';
lambdaOptions.k = 10;
% lambdaOptions.MSE = 4;
lambdaOptions.maxIts = 100;
lambdaOptions.fPlot = true;

mseThreshold = 0.1;
[d, lambdaTrunc, omegaTrunc] = findMinNumLambdaKLE(mseThreshold, a, l, sigma, lambdaOptions);
