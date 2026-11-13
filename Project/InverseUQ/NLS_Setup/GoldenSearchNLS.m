function [theta, P, nlsDetails] = GoldenSearchNLS(b, a, maxIts, thetaRelTol, nTerms, ...
                mode_estimator, covInflationSF, yMeas, controller, fPrintInfo)
    %{
        Runs bisection method to try to find the best starting IC that
        minimizes the NLS cost function
        b = upper bound on search
        a = lower bound on search

        Only searches over constant term because assuming this drives the
        response and dont have any other intuition about other terms.
    %}
       
    % Assumes unimodal function
    phi = 2/(sqrt(5) - 1);
    
    thetaEstPrev = b;
    thetaEst = (b+a)/2;
    itCounter = 0;
    
     while (abs((thetaEst - thetaEstPrev)/thetaEstPrev) > thetaRelTol) ...
             && (itCounter < maxIts)
        c = b - (b-a)/phi;
        d = a + (b-a)/phi;

        if fPrintInfo
            fprintf('Iteration: %i\n', itCounter);
            fprintf('Theta_0 Estimate: %3.3f\n', thetaEst);
        end

        % Run NLS
        [~, ~, nlsDetails_C] = ...
                    RunNLS(mode_estimator, getTheta(c,nTerms), ...
                        covInflationSF, yMeas, controller, false);

        [~, ~, nlsDetails_D] = ...
                    RunNLS(mode_estimator, getTheta(d,nTerms), ...
                        covInflationSF, yMeas, controller,false);

        fC = nlsDetails_C.Jnls(end);
        fD = nlsDetails_D.Jnls(end);


        if fC < fD
            b = d;
        else
            a = c;
        end      

        thetaEstPrev = thetaEst;
        thetaEst = (b+a)/2;
        itCounter = itCounter + 1;    
     end

     % Evaluate final value
     [theta, P, nlsDetails] = ...
                    RunNLS(mode_estimator, getTheta(thetaEst,nTerms), ...
                        covInflationSF, yMeas, controller, false);
end

function theta = getTheta(theta0,nTerms)
    theta = [theta0; zeros(nTerms-1,1)];
end