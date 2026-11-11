function [const] = ModelParams(varargin)
    FEET_TO_M = 0.3048;
    LB_TO_KG = 0.453;

    SLUGS_PerFt3_To_Kg_perM3 = 515.37881839;

    const = struct();

    const.dT = 1E-2;
    const.tEnd = 11; % 10 seconds of ring out time after 1 second step
    const.continuous = true;
    
    const.nStates = 5;
    const.stateNames = {'V','gamma','q','theta','h'};
    const.stateUnits = {'m/s','rad','rad/s','rad','m'};
    const.statePlottingUnits = {'m/s','deg','deg/s','deg','km'};
    const.statePlottingSF = [1, rad2deg(1), rad2deg(1), rad2deg(1), 1E-3];

    const.nInputs = 1;
    const.inputNames = {'deltaE'};
    const.inputUnits = {'rad',};
    const.inputPlottingUnits = {'deg'};
    const.inputPlottingSF = [rad2deg(1)];

    const.g = 9.81; % m/s^2

    const.atmo.fUseSimpleAtmo = false; % Force full atmo in case want to perturb later
    const.atmo.h0 = 8.5E4 * FEET_TO_M; 
    const.atmo.rho0 = 6.7429E-5 * SLUGS_PerFt3_To_Kg_perM3; 
    const.atmo.hs = 2.1358E4 * FEET_TO_M; 

    %{
        NOTE: using "Control Orientated Modelling of an Air-Breathing
        Hypersonic Vehicle" by Parker, Bolender and Doman as the master
        source of truth for simulation parameters. 

        Assuming when any units are normalized by /ft it means they are
        scaled by the vehicle length.

        Assuming that the canard is identical to the elevator in lift/drag.

        Lordy what I'd give for metric units that are consistent between
        papers.
    %}

    vehicleLength = 100*FEET_TO_M; % 100 ft to m
    const.gravimetrics.m = 4377*vehicleLength; % kg (300 lb/ft - assuming normalized by length?)
    const.gravimetrics.Iyy = 5E5*1.35581795*vehicleLength; % kg.m^2 (5E5 slugs.ft^2/ft)
    
    % Aero Coeffs
    const.aero.refArea = (17*FEET_TO_M^2)/FEET_TO_M * vehicleLength; % 17ft^2/ft normalized to vehicle length
    const.aero.refLength = 17*FEET_TO_M; % m 
    
    % NOTE: I bastardised the paper "Control Orientated Modelling of
    % Air-Breathing Hypersonic Vehicles" by Parker et al.
    const.aero.CL_alpha = 4.6773; % rad^-1
    const.aero.CL_deltaE = 7.6224E-1; % rad^-1
    const.aero.CL_0 = -1.8714E-2;

    const.aero.CD_alphaSq = 5.8224; % rad^-2
    const.aero.CD_alpha = -4.5315E-2; % rad^-1
    const.aero.CD_deltaESq = 8.1993E-1; % rad^-2
    const.aero.CD_deltaE = 2.7699E-4; % rad^-1
    const.aero.CD_0 = 1.0131E-2; 

    const.aero.CM_alphaSq = -6.2926; % rad^-2
    const.aero.CM_alpha = -2.1335; % rad^-1
    const.aero.CM_0 = 1.8979E-1;
    const.aero.CM_deltaE = -1.2897; % rad^-1

    %% Dynamic Derivatives
    const.aero.dynDeriv.model = 'off'; %{'off','static','dynamic'}
    % Stabilizing: Took 1/10 of Cmalpha and scaled by 2*Vnom/cbar (~500)
    % Normal value might be closer to -5?
    const.aero.dynDeriv.CM_q_static = -0.2*400; 

    %% Sensor Model
    % TODO for inverse problem

    %% Constraints
    const.constraint.deltaEBounds = [-deg2rad(30), deg2rad(30)];

    %% Controller
    const.controller.sysID.fActive = true;
    const.controller.sysID.mode = 'sine'; % {'step','sine'}

    % Step of fin deflections (zero to max to min to zero) - e.g. impulse
    % response ID)
    const.controller.sysID.stepDuration = 0.1; % s

    % Sine sweep of fin deflections
    const.controller.sysID.sineDuration = 1; % s
    const.controller.sysID.sineFreq = 1; % Hz

    %% Initial Conditions (Default values)
    const.xIC.V =  1500;
    const.xIC.gamma = 0; 
    const.xIC.q = 0;
    const.xIC.theta = 0;
    const.xIC.h = 28000;

    %% Gaussian Process for CM_Alpha
    const.aero.GP.fActive = false;
    const.aero.GP.dataKLE = [];
    const.aero.GP.Y = [];



     %% Pull out specific parameter if required
    nArgs = length(varargin);
    if nArgs > 0
        for iArg = 1:nArgs
            const = const.(varargin{iArg});
        end
    end


end