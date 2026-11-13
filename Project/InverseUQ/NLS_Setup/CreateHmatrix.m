function H = CreateHmatrix(theta, const, controller)

%{
    Assumes we only want to estimate Cm_alpha. 
    Assumes a polynomial model for Cm_alpha and theta are the coefficients
    in increasing powers of time.
%}
    times = 0:const.dT:const.tEnd;
    nTimes = numel(times);
    nMeasPerTime = size(const.measModel.C,1);
    nParamToEst = numel(theta);

    H = zeros(nTimes*nMeasPerTime, nParamToEst);
        
    % integrate up coupled total system
    xIC = [const.xIC.V;
           const.xIC.gamma;
           const.xIC.q;
           const.xIC.theta;
           const.xIC.h];

    S0 = zeros(5,1);
    xAugIC = [xIC; S0];
    const.polyCM_alpha.coeff = theta;
    options = odeset('OutputFcn',@ode45OutputFunc,'RelTol',const.odeOpts.relTol);
    [~, xAugResults] = ode45(@(t,xAug) DynPlusJacobian(t,xAug,controller, const, theta), ...
                            times, xAugIC, options);

    % Extract and build Jacobian
    dCMalpha_dTheta = zeros(nParamToEst,nTimes);
    for i = 1:nParamToEst
        dCMalpha_dTheta(i,:) = times.^(i-1);
    end
    
    for iTime = 1:nTimes
        xAugk = xAugResults(iTime,:)';
        Hk = xAugk(6:end) * dCMalpha_dTheta(:,iTime)';

        tStartIdx = nMeasPerTime*(iTime - 1) + 1;
        tEndIdx = tStartIdx + nMeasPerTime - 1;
        H(tStartIdx:tEndIdx, :) = Hk;
    end
end
