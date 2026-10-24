%{
    Testing to see what i can get it to do
%}

const = ModelParams();

tEnd = 50;
xIC = [1500; 0; 0; 0; 28000]; % Guess

PDControl = @(t,x,const) BasicController(t,x,const);

simData_CBFon = Simulation(xIC, PDControl, tEnd);



plotSimulationResults({simData_CBFon}, {'test'});



%% Controller
function [u, exit] = BasicController(t,x,const)
    
    exit = 0;

    f = 1/10;
    omega = 2*pi*f;
    A = deg2rad(10);
    alphaDes = A * sin(omega*t);
    alphaRateDes = -A*omega*cos(omega*t);

    % Do alpha control
    Kd = 100;
    Kp = 10; % just a guess lol

    alpha = calcAlpha(x);

    % PD Controller
    virtualElevator =  (-Kd*(x(3) - alphaRateDes) - Kp*(alpha - alphaDes));
    deltaE_cmd = min(max(virtualElevator, const.constraint.deltaEBounds(1)), ...
                    const.constraint.deltaEBounds(2));
    
    deltaE_cmd = virtualElevator;
    
    u = [0; deltaE_cmd; 0]; % No engine and no canard
end