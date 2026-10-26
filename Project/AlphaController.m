function [u, exit] = AlphaController(t, x, const, fRefAlpha)
    %{
        Alpha Controller/SysID
    %}
    exit = 0;

    % PD Controller
    if ~const.controller.sysIDActive 
        [alphaDes, alphaRateDes] = fRefAlpha(t);
    
        % Do alpha control
        Kd = 20;
        Kp = 60; % just a guess lol
    
        alpha = calcAlpha(x);
    
        % PD Controller    
        virtualElevator =  (-Kd*(x(3) - alphaRateDes) - Kp*(alpha - alphaDes));

    else
        % Injecting Elevator Doublet
        if t <= 0.1
            virtualElevator = const.constraint.deltaEBounds(1);
        elseif t> 0.1  && t < 0.2
            virtualElevator = const.constraint.deltaEBounds(2);
        else
            virtualElevator = 0;
        end
    end

    % Saturation
    deltaE_cmd = min(max(-virtualElevator, const.constraint.deltaEBounds(1)), ...
                    const.constraint.deltaEBounds(2)); 
       
    % Build Inputs
    u = [0; deltaE_cmd; 0]; % No engine and no canard

end