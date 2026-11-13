function [CL,CD,CM,CMq] = getAeroCoeffs(alpha,deltaE,const, t, x)
    c = const.aero;

    % Handle cm_alpha as a stochastic process - less invasive than a
    % prostate exam but getting close....
    if c.GP.fActive
        %NOTE: might change from time variable to some psuedo heat load
        %variable in future hence passing in state (x) as well as time.
        CM_alpha = getRV_KLE(t, c.GP.dataKLE, c.GP.Y);        
      
    elseif c.movMeanCM_alpha.fActive  % Linear model in alpha
        CM_alpha = c.CM_alpha + c.movMeanCM_alpha.deriv * t;   
       
    elseif c.polyCM_alpha.fActive % Polynomial alpha model
        CM_alpha = 0;
        for i = 1:numel(c.polyCM_alpha.coeff)
            CM_alpha = CM_alpha + c.polyCM_alpha.coeff(i)*t.^(i-1);
        end  

    else % Default case
        CM_alpha = c.CM_alpha;
    end

    % Standard Static Coefficients vs Time
    CL = c.CL_alpha*alpha + c.CL_deltaE*deltaE + c.CL_0;

    CD = c.CD_alphaSq*alpha.^2 + c.CD_alpha*alpha + ...
        c.CD_deltaE*deltaE + c.CD_deltaESq*deltaE.^2 + ...
        c.CD_0;

    CM = c.CM_alphaSq*alpha.^2 + CM_alpha*alpha + c.CM_deltaE*deltaE + ...
         c.CM_0;

    % Dynamic derivatives
    if strcmpi('off',c.dynDeriv.model)
        CMq = 0;
    elseif strcmpi('static',c.dynDeriv.model)
        CMq = c.dynDeriv.CM_q_static;
    else
        error('invalid dynamic derivative option');
    end    
end