function [simData] = Simulation(fController, const)
    %{
        fController = function handler to a controller of the form:
            u = f(t,x,const), where x = state of vehicle (V, gamma, q, theta,
            height), u = [deltaE], const=ModelParams(). Using time as an option too
            in case trying to track a reference trajectory embedded in
            fController.
    %}
    xIC = [const.xIC.V;
           const.xIC.gamma;
           const.xIC.q;
           const.xIC.theta;
           const.xIC.h];

    %% Data Packaging
    simData = struct();
    simData.times = 0:const.dT:const.tEnd;
    nTimes = numel(simData.times);


    %% Run Sim
    xResults = zeros(nTimes, const.nStates);
    if const.continuous == true
        options = odeset('OutputFcn',@ode45OutputFunc);
        [~, xResults] = ode45(@(t,x) plantDynCL(t,x,fController, const), ...
                                simData.times, xIC, options);
    else
        xResults(1,:) = xIC';
        xk = xIC;
        u = zeros(1,1); % Default if controller fails
        for k=2:nTimes
            tkprev = simData.times(k-1);
            tk = simData.times(k);
            [u,~] = fController(tkprev, xk, const);
            [~, xtraj] = ode45(@(t,x) getDynamics(x, u, const), ...
                                    [0,const.dT], xk);
            xk = xtraj(end,:)';
            disp(tk);
            xResults(k,:) = xk';
        end
    end
    simData.x = xResults';   

    %% Extract Control Actions and other Useful Quantities
    simData.u = NaN(const.nInputs,nTimes);
    simData.slackVar = NaN(nTimes,1);
    simData.qBar = NaN(nTimes,1);
    for i = 1:nTimes
        %Control Input
        [simData.u(:,i)] = ...
                    fController(simData.times(i), ...
                                simData.x(:,i), const);

        %Qbar
        rho = getAtmo(simData.x(end,i), const);
        simData.qBar(i) = 0.5*rho*simData.x(1,i)^2;
    end
end

%% Helpers
function xDot = plantDynCL(t, x, fController, const)
    %{
        Wrapper for CL dynamics to make everything smooth and continous.
        Aka "we have ZOH at home"
    %}    
    [u,~] = fController(t, x, const);
    xDot = getDynamics(x, u, const, t);
end