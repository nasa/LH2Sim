% Defines LH2 model parameters
%
% Copyright (c) 2017 United States Government as represented by
%     the Administrator of the National Aeronautics and Space Administration.
%     All Rights Reserved.
%
clear LH2Model;

psiToPa = 6894.75729;

% environmental parameters
LH2Model.pAtm = 1.01325e5;			% [Pa] pressure of atmosphere
LH2Model.QdotEW = 20e3;				% [W] external heat transfer

% Tank parameters
LH2Model.R = 8.41/2;				% [m] radius of external tank
LH2Model.A = pi*(LH2Model.R)^2;		% [m^2] cross section area of external tank
LH2Model.H = 29.5656;				% [m] height of external tank
LH2Model.g = 9.8;					% [m/s^2] acceleration due to gravity
LH2Model.VTotal = LH2Model.A*LH2Model.H;	% [m^3] total volume of ET
LH2Model.cw = 5e2;					% [J/kg/K] thermal inertia term for tank wall
LH2Model.tw = 0.1;					% [m] thickness of wall
LH2Model.rhow = 2700;				% [kg/m^3] density of wall

% LH2 parameters
LH2Model.T_c = 33.2;				% [K] critical temperature of hydrogen
LH2Model.p_c = 1.315e6;				% [Pa] critical pressure of hydrogen
LH2Model.lambda = 5;				% dimensionless exponent for saturated vapor
LH2Model.rhoL = 71.1;				% [kg/m^3] liquid density
LH2Model.c_L = 9450;				% [J/kg/K] specific heat of liquid
LH2Model.kappaL = 0.0984;			% [W/mK] thermal conductivity of liquid
LH2Model.qh0 = 4.47e5;				% [J/kg] specific heat of evaporation at Tc
LH2Model.mu = 3.4e-6;				% [Pa*s] dynamic viscosity vapor
LH2Model.muL = 1.3e-5; 0.01;				% [Pa*s] dynamic viscosity liquid
LH2Model.betaL = 0.02;				% [1/K] coefficient thermal expansion liquid

% GH2 (vapor) parameters
LH2Model.R_v = 4120; 4124;				% [J/kg/K] vapor constant
LH2Model.c_v = 10160; 6490;				% [J/kg/K] specific heat of hydrogen at V=const; rotational degrees are frozen
LH2Model.c_p = 14320; LH2Model.c_v + LH2Model.R_v;
LH2Model.gamma = 5/3;				% dimensionless parameter
LH2Model.Gamma = ((LH2Model.gamma+1)/2)^((LH2Model.gamma+1)/2/(LH2Model.gamma-1)); % dimensionless parameter
LH2Model.kappa_v = 0.0166;			% [W/mK] thermal conductivity of saturated vapor at T=20K

% GHe (gas) parameters
LH2Model.R_g = 2077;				% [J/kg/K] helium gas constant
LH2Model.kappa_g = 0.0262;			% [W/mK] thermal conductivity of helium at T=20K
LH2Model.c_gv = 3121;				% [J/kg/K] specific heat of helium at V=const
LH2Model.c_gp = 5193;				% [J/kg/K] specific heat of helium at P=const

% Stratification parameters
LH2Model.n = 10;					% [] number of vapor/gas layers
LH2Model.k = 10;					% [] number of liquid layers

% Optimal Grid parameters
LH2Model.nL = 3;					% [] grid size (liquid)
LH2Model.nV = 3;					% [] grid size (vapor)
LH2Model.tmin = 0.1;				% [s] grid time constant (vapor)

% Initial conditons
LH2Model.p0 = 2*14.7*psiToPa;
LH2Model.TL0 = 20;
LH2Model.Tv0 = 20;
LH2Model.Tw0 = 21;
LH2Model.TwL0 = LH2Model.Tw0;
LH2Model.Twv0 = LH2Model.Tw0;
LH2Model.mL0 = 6e4;
LH2Model.mv0 = LH2Model.p0/2*(LH2Model.VTotal-LH2Model.mL0/LH2Model.rhoL)/LH2Model.Tv0/LH2Model.R_v;
LH2Model.mg0 = LH2Model.p0/2*(LH2Model.VTotal-LH2Model.mL0/LH2Model.rhoL)/LH2Model.Tv0/LH2Model.R_g;

% Input conditions
LH2Model.JLOut = 10;

% Vent valve parameters
LH2Model.Venting = 0;				% whether venting is enabled or not
LH2Model.VentState = 0;				% state of vent valve
LH2Model.SVent = 0.0005;			% [m^2] vent valve orifice area
%LH2Model.KVent = 1e-6;					% dimensionless valve coefficient
LH2Model.pVentLow = 2*LH2Model.pAtm;%38.7*psiToPa;	% [Pa] vent valve lower pressure threshold
LH2Model.pVentHigh = 2.02*LH2Model.pAtm;%41.7*psiToPa;	% [Pa] vent valve upper pressure threshold
LH2Model.Pressurization = 0;		% whether pressurizing or not (using 3 atm for external pressure)

% Temperature sensor locations
LH2Model.sensorPositions = [2 12 16 28];

% Solver options
LH2Model.tFinal = 1000;				% [s] final time
LH2Model.relTol = 1e-4;				% relative tolerance for ode solver
