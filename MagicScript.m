%%
%%  Author: Tiago Monteiro
%%
%%
close all
clear
clc

addpath('src/');

YES = 1;
NO  = 0;
ON  = 1;
OFF = 0;
Sim.horn_number = 0;
Sim.adapt_number = 0;
%------------------------------------------------------------------------------%
%                                   MagicScript                                %
%------------------------------------------------------------------------------%
%%
%% OpenEMS based app for corrugated Horn antennas simulation. 
%% Edit "user editable parameters" for change horn characteristic. This program
%% can sweep parameters if they are defined as a linear matrix.
%%
%%    E.g. Sim.ao = 5:1:10;    % will sweep between 6 simulations with ao 
%%                                  starting at 5 and finishing at 10.

%
% >>>>>> User Editable parameters >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%

%___ General parameters __________________
RUN_SIMULATION  = ON;
Sim.MAKE_HORN   = OFF;
Sim.MAKE_ADAPT  = ON;

Sim.output_path = 'outputs/';
Sim.Sim_Path    = 'tmp/';
Sim.Sim_CSX     = 'Simulation.xml';

Sim.fmin    = 8;                     % Minimum frequency in GHz
Sim.fmax    = 12;                    % Maximum frequency in GHz
Sim.fcalc   = 10;                    % Frequency to calculate fields

%___ Horn parameters _____________________
Sim.horn_ai      = 22.86;
Sim.horn_bi      = 10.16;

Sim.horn_ao      = 145;%63.75;
Sim.horn_bo      = 120;%49.78;

%__ Corrugated profile ____
%
%               |width| slot |     |  corr_step |
%                _____        _____        _____  _
%               |     |      |     |      |     |   depth
%           ____|     |______|     |______|     | _
%                                               |
%           ____________________________________|
%
%               slot   = corr_step * delta
%               width  = corr_step * (1-delta)
%
% If Sim.depth is set to 0 it goes from fcalc/2 to fcalc/4 across flare length.
% If Sim.USE_CORRUGATIONS_A/B are set to NO or OFF then depth_a/b are ignored.
%
Sim.horn_corr_step           = 3;
Sim.horn_delta               = 0.75;
% A walls setup
    Sim.HORN_USE_CORRUGATIONS_A  = OFF;
    Sim.horn_depth_a             = 0;
    Sim.horn_a_jump              = 0;        % First corrugation jump
    Sim.HORN_PROFILE_FOR_A       = 1;        % 1=Linear, 2=Tangential, 3=Exponential, 4=Two phased linear.
    % Only for Tangential profile
    Sim.horn_A_wall_tan_A        = 1;        % A coeficient for Tangential wall A
    Sim.horn_A_wall_tan_rho      = 2;        % rho coeficient for Tangential wall A
    % Only fot two phase linear
    Sim.horn_first_a_length      = 0.1;      % First length in horn flare
    Sim.horn_first_ao            = 70;       % First flare aperture for a walls

% B walls setup
    Sim.HORN_USE_CORRUGATIONS_B  = ON;
    Sim.horn_depth_b             = 0;
    Sim.horn_b_jump              = 2;        % First corrugation jump
    Sim.HORN_PROFILE_FOR_B       = 4;        % 1=Linear, 2=Tangential, 3=Exponential, 4=Two phased linear.
    % Only for Tangential profile
    Sim.horn_B_wall_tan_A        = 1;        % A coeficient for Tangential wall B
    Sim.horn_B_wall_tan_rho      = 2;        % rho coeficient for Tangential wall B
    % Only fot two phase linear
    Sim.horn_first_b_length      = 0.2;      % First length in horn flare
    Sim.horn_first_bo            = 25;       % First flare aperture for b walls

Sim.horn_wg_length           = 60;           % Length of feeding waveguide
Sim.horn_num_of_corrugations = 35;           % Number of corrugations
Sim.horn_straight_width      = 2;
Sim.horn_cap_width           = 2;

Sim.horn_exc_mode            = 'TE10';       % Port excitation mode

Sim.HORN_SHOW_STRUCTURE_FIGURES      = YES;
Sim.HORN_SUBSTRACT_LEFTOVERS         = YES;  % Subtracts horn leftovers with air volume

Sim.HORN_TIME_STEPS = 10000;
Sim.horn_n_cell     = 40;                    % cell_size = lambda_fmax / n_cell


%___ Adapter parameters __________________
Sim.ADAPT_TIME_STEPS = 50000; %max. number of timesteps

% waveguide dimensions and mode
Sim.adapt_m = 1;
Sim.adapt_n = 0;

Sim.adapt_length          = 60;
Sim.adapt_b               = 10.16;            % Waveguide width 
Sim.adapt_a               = 22.86;            % Waveguide heigth
Sim.adapt_WallThickness   = 2;                % Walls thickness
Sim.adapt_BackShort       = 7;                % Distance from short to center of probe


%%%%% Conector N %%%%%
Sim.adapt_InnerCond_N           = 3.15;             % Inner diameter
Sim.adapt_OuterCond_N           = 8.13;             % Inner diam of outer conductor
Sim.adapt_OuterCondOD_N         = 15.8;             % Outer diam of outer conductor
Sim.adapt_ProbeDepth            = 5.5;              % Probe insertion depth inside waveguide
Sim.adapt_ProbeRad              = 1;                % Probe rad
Sim.adapt_dielectric_intrusion  = 2;             % Dielectric intrusion on waveguide
Sim.ADAPT_ADD_SPHERE            = NO;                   % Add probe end sphere
Sim.adapt_sph_rad               = Sim.adapt_ProbeRad;   % Probe end sphere radius
Sim.adapt_N_Length              = 10.72;            % Length of N connector
Sim.adapt_epsR                  = 2.08;             % Teflon permitivity

Sim.adapt_mesh_res        = [.3 .3 .3];
Sim.adapt_space           = 5;
Sim.adapt_exc_mode        = 'TE10';

%
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< End of User Editable parameters <<<<<<
%

simulation_core(Sim, RUN_SIMULATION);

disp(">>-------------------- Simulation fineshed! --------------------<<");
