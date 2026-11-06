function [CL,CD,CM,CMq] = getAeroCoeffs(alpha,deltaE,const)
    c = const.aero;

    CL = c.CL_alpha*alpha + c.CL_deltaE*deltaE + c.CL_0;

    CD = c.CD_alphaSq*alpha.^2 + c.CD_alpha*alpha + ...
        c.CD_deltaE*deltaE + c.CD_deltaESq*deltaE.^2 + ...
        c.CD_0;

    CM = c.CM_alphaSq*alpha.^2 + c.CM_alpha*alpha + c.CM_deltaE*deltaE + ...
         c.CM_0;

    % Dynamic derivatives
    if strcmpi('off',c.dynDeriv.model)
        CMq = 0;
    elseif strcmpi('static',c.dynDeriv.model)
        CMq = c.dynDeriv.CM_q_static;
    else
        error('invalid dynamic derivative option');
    end
    
    % Make statically stable
    if c.fMakeStaticallyStable
        CM = -c.CM_alphaSq*alpha.^2 + -c.CM_alpha*alpha + c.CM_deltaE*deltaE + ...
            c.CM_0;
    end
    
end