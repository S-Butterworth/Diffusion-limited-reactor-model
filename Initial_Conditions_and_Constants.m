function [C,F,data] = Initial_Conditions_and_Constants()

%% Loading data for eta based on temperature and days on stream
iptDir = fileparts(mfilename('fullpath'));
outDirdata = fullfile(iptDir, 'data');
if ~isfolder(outDirdata)
    mkdir(outDirdata);
end
data = load(fullfile(outDirdata,'pellet_eta_table_activity.mat'));

% Product requiring annual production rates 
selec_prod = 3;
stoich = [-1,-1/2,1,0];
inert_ratio = 0.79/0.21;
inert_ref = 2;

mr = [28.053/1000,31.999/1000,44.0526/1000,28/1000];

OPDays = 350;      % Operational Days
C.NoReactors = 20; % Number of Reactors
C.NoTubes = 100;   % Number of Tubes per Reactor
Seconds_in_day = 24*60^2;

% C2H4O kgyear-1 (350 days operational time), converted to kgs-1 per reactor
Totoal_Prod = 160e6; 
Prod = Totoal_Prod / (OPDays*Seconds_in_day*C.NoReactors);

C.X = 0.75; % Conversion rate
Fexit = Prod / mr(selec_prod);

% Determines initial rate of flow

C.F = (sum(-1*stoich,1).*(sum(-1*stoich,1)>0)).*Fexit ./ C.X;

inert = inert_ratio*sum(-1*stoich(inert_ref),1);
C.F(end) = inert*Fexit / C.X;

MF = C.F.*mr;
C.x = C.F/sum(C.F);

PT = 1E6;

% NASA polynomials A,B,C,I
NASA = [[ 2.6E+00, 3.2E+00, 3.8E+00, 3.5E+00]; ... % a0
    [ 1.3E-02, 1.7E-03,-9.4E-03,-1.2E-04]; ... % a1
    [-5.5E-06,-7.7E-07, 8.0E-05,-5.0E-07]; ... % a2
    [ 1.1E-09, 1.5E-10,-1.0E-07, 2.4E-09]; ... % a3
    [-7.5E-14,-1.1E-14, 4.0E-11,-1.4E-12]; ... % a4
    [ 4.8E+03,-1.0E+03,-7.6E+03,-1.0E+03]; ... % a5
    [ 7.6E+00, 6.2E+00, 7.8E+00, 3.0E+00]];    % a6

%% Constants
C.k =  2.54E-08;                            % Pre-exponential (reference) rate constant, evaluated at T0
C.F0 = sum(C.F);                            % Initial total molar flow rate (sum across all species)
C.PT = PT;                                  % Total system pressure, Pa
C.Td = 35E-3;                               % Reactor tube diameter, m
C.CSA = pi() * (C.Td/2)^2 * C.NoTubes ;     % Reactor cross-sectional area
C.vf = 0.55;                                % Bed voidage (void fraction)
C.pc = 1650;                                % Catalyst bulk density, kg/m^3
C.G = sum(MF)/C.CSA;                        % Mass flux (superficial mass velocity), kg/(m^2 s)
C.dp = 5.0E-3;                              % Catalyst pellet diameter, m
C.dens = 6.616;                             % Fluid density, kg/m^3
C.mu = 2.782E-5;                            % Fluid viscosity, Pa.s
C.NASA = NASA;                              % NASA polynomial coefficients (for Cp/enthalpy evaluation)
C.EA = 39E3;                                % Activation energy, J/mol
C.R = 8.31446;                              % Universal gas constant, J/(mol K)
C.T0 = 200;                                 % Inlet/reference temperature, deg C
C.A = C.k*exp(C.EA/(C.R*(C.T0+273.15)));    % Arrhenius pre-exponential factor, referenced to T0
C.a0 = 1;                                   % Initial (fresh) catalyst activity, normalized
C.kd = 5.20E-2;                             % Deactivation rate constant
C.t = 0;                                    % Catalyst time-in-use, days (0 = fresh catalyst)
C.n = 2;                                    % Deactivation order
C.stoich = stoich;                          % Reaction stoichiometric coefficients
C.coeff = [1/3,2/3,0,0];                    % Reaction order (exponent) for each species in the rate law
C.D = [1.556E-07, 2.831E-06];               % Effective diffusivities for species A and B, m^2/s
C.cat = 100;                                % Initial catalyst mass/loading

C.threshold = 1e-6;          % Dead-core detection threshold on C_A (triggers switch attempt)
C.minRcToSwitch = 0.01;      % Minimum estimated r_c required before attempting the free-boundary switch
C.RelTol1 = 1e-5;            % bvp4c relative tolerance — Stage 1 (regularized, full-domain) solve
C.RelTol2 = 1e-3;            % bvp4c relative tolerance — Stage 2 (free-boundary) solve

%% Flags
F.thermo = true;    % True - adiabatic, False - isothermal
F.activity = true;  % True, false
F.diffusion = true; % True, false
F.pressure = true;  % True, false - isobaric
F.stats = true;     % Set up information
end