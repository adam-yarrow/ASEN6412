%{
    TODO - README
%}

%% Testing
% Check root finding
omega = 0:0.001:20;
a = 0.5;
l = 2;

omegaSol_odd = [];
omegaSol_Even = [];
omegaSolTanVal_odd = [];
omegaSolTanVal_even = [];
for i = 1:10
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
hold on;
plot(omega, tan(omega*a),'k');
plot(omega, 1./(l*omega), 'r--');
plot(omega, -l*omega,'b--');
plot(omegaSol_Even, omegaSolTanVal_even, 'bo');
plot(omegaSol_odd, omegaSolTanVal_odd, 'ro');
legend('tan','odd','even','solutions');
ylim([-60,20]);
grid on;

