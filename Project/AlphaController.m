function [u, exit] = AlphaController(t, x, const, fRefAlpha)
    %{
        PD Controller to 
    %}
    exit = 0;

    % f = 1/10;
    % omega = 2*pi*f;
    % A = deg2rad(10);
    % alphaDes = A * sin(omega*t);
    % alphaRateDes = -A*omega*cos(omega*t);

    [alphaDes, alphaRateDes] = fRefAlpha(t);

    alphaDes = 0;
    alphaRateDes = 0;

    % Do alpha control
    Kd = 20;
    Kp = 60; % just a guess lol

    alpha = calcAlpha(x);

    % PD Controller
    virtualElevator =  (-Kd*(x(3) - alphaRateDes) - Kp*(alpha - alphaDes));

    % Saturation
    deltaE_cmd = min(max(-virtualElevator, const.constraint.deltaEBounds(1)), ...
                    const.constraint.deltaEBounds(2)); 
        
    % Build Inputs
    u = [0; deltaE_cmd; 0]; % No engine and no canard

end