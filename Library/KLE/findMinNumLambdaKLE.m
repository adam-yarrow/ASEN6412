function [d, lambda, omega] = findMinNumLambdaKLE(mseThreshold, a, l, sigma, options)
%{
    Calculates the minimum number of eigenfunctions required to obtain an
    approximation of a random process which achieves a specified Mean
    Square Error threshold. 

    Assumes an exponentially decaying covariance function.

    Will return the minimum number of lambda and their values (+omega
    values);

    mseThreshold = threshold as a fraction for acceptable error
    a = domain bound
    l = characteristic length of exp cov function
    sigma = std dev for exp cov function
    options = struct containing
        - mode = 'analytical','approximate'
        - k = (Optional) If approximate, k is number of modes in advance to estimate
        with
        - MSE = (Optional) If analytical this is the analytical MSE.
%}

alpha = 1 - mseThreshold;
mode = options.mode;


% Analytical Approach
if strcmpi('analytical',mode)
    [d, lambda, omega] = analyticalSol(alpha, a, l, sigma, options);
elseif strcmpi('approximate',mode)
    [d, lambda, omega] = approxSol(alpha, a, l, sigma, options);
else
    funCallStack = dbstack;
    methodName = funCallStack(1).name;
    error('Invalid mode provided for function: %s', methodName);
end   

end


%% Helpers
function [d, lambda, omega] = approxSol(alpha, a, l, sigma, options)
truncErrorFrac = 0;
lambda = [];
omega = [];
i = 1;

lambdaSum = 0;
lambdaSumExtended = 0;

while truncErrorFrac < alpha
        
    if i == 1
        for j = 1:options.k + 1
            omega_j = getOmegaKLE(j,a,l);
            lambda_j = getLambdaKLE(l,omega_j,sigma);
            lambda = [lambda, lambda_j];
            omega = [omega, omega_j];

            lambdaSumExtended = lambdaSumExtended + lambda_j;
            if j == 1
                lambdaSum = lambdaSum + lambda_j;
            end
        end
    else
        % only add the next in the series (e.g. i + k)
        omega_j = getOmegaKLE(i+options.k,a,l);
        lambda_j = getLambdaKLE(l,omega_j,sigma);

        lambda = [lambda, lambda_j];
        omega = [omega, omega_j];

        lambdaSumExtended = lambdaSumExtended + lambda_j;
        lambdaSum = lambdaSum + lambda(i);
    end

    truncErrorFrac = lambdaSum / lambdaSumExtended;
    i = i + 1;
end

d = i - 1;

% Plotting
if options.fPlot
    figure('Name','Eigenvalues for Exponential Covariance Function');
    plot(1:numel(lambda), lambda ./ sum(lambda),'r*','MarkerSize',8);
    grid on;
    xlabel('i','FontSize',14)
    txt = '$\lambda_i $ / $ \sum_{j=1}^{k}\lambda_j$';
    ylabel(txt,'Interpreter','latex','FontSize',14);
    title('Normalized Eigenvalues vs Mode Number','FontSize',14)
end
end

function [d, lambda, omega] = analyticalSol(alpha, a, l, sigma, options)
truncErrorFrac = 0;
i = 1;

lambda = [];
omega = [];

lambdaSum = 0;
% Iterate until meet threshold
while (truncErrorFrac < alpha)
    omega_i = getOmegaKLE(i,a,l);
    lambda_i = getLambdaKLE(l,omega_i,sigma);
    lambdaSum = lambdaSum + lambda_i;

    truncErrorFrac = lambdaSum / options.MSE;
       
    lambda = [lambda, lambda_i];
    omega = [omega, omega_i];
    
    i = i+1;

    % Break if exceed max its
    if i > options.maxIts
        error('Could not converge to valid minimum d')
    end
end

d = i -1;

% Plotting
if options.fPlot
    figure('Name','Eigenvalues for Exponential Covariance Function');
    plot(1:numel(lambda), lambda ./ options.MSE,'r*','MarkerSize',8);
    grid on;
    xlabel('i','FontSize',14)
    txt = '$\lambda_i $ / $ <\int X(t,\omega) dt> $';
    ylabel(txt,'Interpreter','latex','FontSize',14);
    title('Normalized Eigenvalues vs Mode Number','FontSize',14)
end
end


