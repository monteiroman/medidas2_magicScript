%%
%%  Author: Tiago Monteiro
%%
%% Clean workspace
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
    Sim.horn_a_jump              = 0;
    Sim.HORN_PROFILE_FOR_A       = 1;        % 1=Linear, 2=Tangential, 3=Exponential, 4=Two phased linear.
    % Only for Tangential profile
    Sim.horn_A_wall_tan_A        = 1;        % A coeficient for Tangential wall A
    Sim.horn_A_wall_tan_rho      = 2;        % rho coeficient for Tangential wall A
    % Only fot two phase linear
    Sim.horn_first_a_length      = 0.1;
    Sim.horn_first_ao            = 70;

% B walls setup
    Sim.HORN_USE_CORRUGATIONS_B  = ON;
    Sim.horn_depth_b             = 0;
    Sim.horn_b_jump              = 3.5;
    Sim.HORN_PROFILE_FOR_B       = 4;        % 1=Linear, 2=Tangential, 3=Exponential, 4=Two phased linear.
    % Only for Tangential profile
    Sim.horn_B_wall_tan_A        = 1;        % A coeficient for Tangential wall B
    Sim.horn_B_wall_tan_rho      = 2;        % rho coeficient for Tangential wall B
    % Only fot two phase linear
    Sim.horn_first_b_length      = 0.2;
    Sim.horn_first_bo            = 25;

Sim.horn_wg_length           = 60;           % Length of feeding waveguide
Sim.horn_num_of_corrugations = 35;
Sim.horn_straight_width      = 2;
Sim.horn_cap_width           = 2;

Sim.horn_exc_mode            = 'TE10';       % Port excitation mode

Sim.HORN_SHOW_STRUCTURE_FIGURES      = YES;
Sim.HORN_SUBSTRACT_LEFTOVERS         = YES;  % Subtracts horn leftovers with air volume

Sim.HORN_TIME_STEPS = 10000;
Sim.horn_n_cell     = 40;                   % cell_size = lambda_fmax / n_cell


%___ Adapter parameters __________________
Sim.ADAPT_TIME_STEPS = 50000; %max. number of timesteps

% waveguide dimensions and mode
Sim.adapt_m = 1;
Sim.adapt_n = 0;

Sim.adapt_length          = 60;
Sim.adapt_b               = 10.16;            % waveguide width 
Sim.adapt_a               = 22.86;            % waveguide heigth
Sim.adapt_WallThickness   = 2;                % walls thickness
Sim.adapt_BackShort       = 7;                % distance from short to center of probe


%%%%% Conector N %%%%%
Sim.adapt_InnerCond_N           = 3.04;            %inner diameter
Sim.adapt_OuterCond_N           = 8.13;            %inner diam of outer conductor
Sim.adapt_OuterCondOD_N         = 15.8;               % outer diam of outer conductor
Sim.adapt_ProbeDepth            = 5;              % Probe insertion depth inside waveguide
Sim.adapt_dielectric_intrusion  = 1:0.5:4;
Sim.ADAPT_ADD_SPHERE            = NO;
Sim.adapt_sph_rad               = 2;
Sim.adapt_N_Length              = 10.72;            % length of N connector
Sim.adapt_epsR                  = 2.08;             % Teflon permitivity

Sim.adapt_mesh_res        = [.5 .5 .5];
Sim.adapt_space           = 5;
Sim.adapt_exc_mode        = 'TE10';

%
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< End of User Editable parameters <<<<<<
%

simulation_core(Sim, RUN_SIMULATION);

disp(">>-------------------- Simulation fineshed! --------------------<<");
