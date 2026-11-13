function plotNLS(nlsDetails,xTrue)
% Diagnostics Plots
figure('Name','Diagnostics Plot');

subplot(5,1,1);
relJnls = abs((nlsDetails.Jnls(2:end) - nlsDetails.Jnls(1:end-1)) ...
    ./ nlsDetails.Jnls(1:end-1));
semilogy(2:nlsDetails.nIts, relJnls);
xlabel('Iteration');
ylabel('Jnls Relative Reduction');
grid on;

subplot(5,1,2);
plot(nlsDetails.alpha);
xlabel('Iteration');
ylabel('Alpha [0,1]');
grid on;

subplot(5,1,3);
plot(nlsDetails.nAlphaItPerIt);
xlabel('Iteration');
ylabel('Alpha Iterations');
grid on;

subplot(5,1,4);
semilogy((sum((nlsDetails.x0Log - xTrue(:,1)).^2,1).^0.5));
xlabel('Iteration');
ylabel('||x_{k} - x_{true}||_2')
grid on;

subplot(5,1,5);
x0Convergence = sum((nlsDetails.x0Log(:,2:end) - nlsDetails.x0Log(:,1:end-1)).^2,1).^0.5;
semilogy(2:1:numel(x0Convergence)+1, x0Convergence);
xlabel('Iteration');
ylabel('||x_{k} - x_{k-1}||_2');
grid on;

sgtitle('Diagnostics Information');

end