function [u, exit] = AlphaController(t, x, const, fRefAlpha)
    %{
        Alpha Controller/SysID
    %}
    exit = 0;

    % PD Controller
    if ~const.controller.sysID.fActive 
        [alphaDes, alphaRateDes] = fRefAlpha(t);
    
        % Do alpha control
        Kd = 20;
        Kp = 60; % just a guess lol
    
        alpha = calcAlpha(x);
    
        % PD Controller    
        virtualElevator =  (-Kd*(x(3) - alphaRateDes) - Kp*(alpha - alphaDes));

    else
        sysId = const.controller.sysID;
        % Injecting Elevator Doublet
        if t <= sysId.duration
            virtualElevator = const.constraint.deltaEBounds(1);
        elseif t > sysId.duration  && t < 2*sysId.duration
            virtualElevator = const.constraint.deltaEBounds(2);
        else
            virtualElevator = 0;
        end
    end

    % Saturation
    deltaE_cmd = min(max(-virtualElevator, const.constraint.deltaEBounds(1)), ...
                    const.constraint.deltaEBounds(2)); 
       
    % Build Inputs
    u = deltaE_cmd;

end