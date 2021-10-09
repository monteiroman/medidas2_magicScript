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

%%
%%  Simulation function
%%
function simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION)
    %%  
    %%  This function makes the simulation directory, calls make_horn for make 
    %%  the 3D structure and then calls run_simulation if it is required to. 
    %%  
    %%  Parameters:
    %%          Sim_Path: Generic simulation path for the simulation data.
    %%          Sim_CSX: OpenEMS XML file name.
    %%          Sim: Horn simulation parameters.
    %%          sweep_type: String detailing sweep type.
    %%          RUN_SIMULATION: Flag that determines if the simulation must be 
    %%                          done.
    %%  Returns:
    %%          void
    %%

    %% ----->> Output simulation particular folder <<-----
    sweep_type = strcat('/', sweep_type);
    if (Sim.horn_number == 0);
        sweep_path = strcat(sweep_type, '_Horn_output');
    else
        sweep_path = strcat(sweep_type, sprintf('_Horn_n%d_output',Sim.horn_number));
    endif
    Sim.output_path = strcat(Sim.output_path, sweep_path);
    [status, message, messageid] = mkdir(Sim.output_path);

    %% ----->> Make horn structure and save openEMS XML simulation file at "<Sim_Path>/<Sim_CSX>"  <<-----
    [port, nf2ff] = make_horn(Sim_Path, Sim_CSX, Sim);

    %% ----->> Simulate and store 3D files <<-----
    if (RUN_SIMULATION);
        run_simulation(Sim_Path, Sim_CSX, Sim, port, nf2ff);
        movefile(strcat(Sim_Path, '/*.stl'), Sim.output_path);
        movefile(strcat(Sim_Path, '/*.vtk'), Sim.output_path);
        movefile(strcat(Sim_Path, '/*.txt'), Sim.output_path);
    endif
endfunction

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

% ______ User Editable parameters ______________________________________________
RUN_SIMULATION = ON;

output_path = 'outputs/';
Sim_Path    = 'tmp/';
Sim_CSX     = 'Corrugated_Horn.xml';

Sim.fmin    = 8;                     % Minimum frequency in GHz
Sim.fmax    = 12;                    % Maximum frequency in GHz
Sim.fcalc   = 10;                    % Frequency to calculate fields

Sim.ai      = 22.86;
Sim.bi      = 10.16;

Sim.ao      = 145;%63.75;
Sim.bo      = 120;%49.78;

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
% If Sim.USE_CORRUGATIONS_A/B are set to NO or OFFthen depth_a/b are ignored.
%
Sim.corr_step           = 6;
Sim.delta               = 0.65:0.05:0.75;
Sim.USE_CORRUGATIONS_A  = OFF;
Sim.depth_a             = 0;
Sim.USE_CORRUGATIONS_B  = ON;
Sim.depth_b             = 0;

Sim.wg_length           = 60;           % Length of feeding waveguide
Sim.num_of_corrugations = 35;
Sim.straight_width      = 2;
Sim.cap_width           = 2;

Sim.exc_mode            = 'TE10';

Sim.SHOW_STRUCTURE_FIGURES      = YES;
Sim.SUBSTRACT_LEFTOVERS         = YES;

Sim.TIME_STEPS  = 10000;
Sim.n_cell      = 20;                   % cell size: lambda/n_cell

Sim.USE_PROFILE = 1;                    % 1=Linear, 2=Tangential, 3=Exponential
Sim.output_path = output_path;

% ______ End of User Editable parameters _______________________________________


%% ----->> Generic simulation output folder <<----- 
confirm_recursive_rmdir(0);                                 % Do not asks if remove directory
[status, message, messageid] = rmdir( output_path, 's');    % Clear previous directory
[status, message, messageid] = mkdir( output_path );        % Create empty simulation folder

%% ----->> Check parameters and sweep if necessary <<----- 
ai_len = length(Sim.ai);
bi_len = length(Sim.bi);
ao_len = length(Sim.ao);
bo_len = length(Sim.bo);
corr_step_len   = length(Sim.corr_step);
delta_len       = length(Sim.delta);
depth_a_len     = length(Sim.depth_a);
depth_b_len     = length(Sim.depth_b);
wg_length_len   = length(Sim.wg_length);
num_of_corrugations_len = length(Sim.num_of_corrugations);

if (ai_len > 1);
    ai_values = Sim.ai;

    for i = 1:ai_len;
        close all                                   % Prevent memory leakage
        Sim.horn_number = i;
        Sim.ai          = ai_values(i);
        sweep_type      = sprintf('ai_sweep_%.2f', ai_values(i));
        simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
    endfor
elseif (bi_len > 1);
    bi_values = Sim.bi;

    for i = 1:bi_len;
        close all
        Sim.horn_number = i;
        Sim.bi          = bi_values(i);
        sweep_type      = sprintf('bi_sweep_%.2f', bi_values(i));
        simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
    endfor
elseif (ao_len > 1);
    ao_values = Sim.ao;

    for i = 1:ao_len;
        close all
        Sim.horn_number = i;
        Sim.ao          = ao_values(i);
        sweep_type      = sprintf('ao_sweep_%.2f', ao_values(i));
        simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
    endfor
elseif (bo_len > 1);
    bo_values = Sim.bo;

    for i = 1:bo_len;
        close all
        Sim.horn_number = i;
        Sim.bo          = bo_values(i);
        sweep_type      = sprintf('bo_sweep_%.2f', bo_values(i));
        simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
    endfor
elseif (corr_step_len > 1);
    corr_step_values = Sim.corr_step;

    for i = 1:corr_step_len;
        close all
        Sim.horn_number = i;
        Sim.corr_step   = corr_step_values(i);
        sweep_type      = sprintf('corr_step_sweep_%.2f', corr_step_values(i));
        simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
    endfor
elseif (delta_len > 1);
    delta_values = Sim.delta;

    for i = 1:delta_len;
        close all
        Sim.horn_number = i;
        Sim.delta       = delta_values(i);
        sweep_type      = sprintf('delta_sweep_%.2f', delta_values(i));
        simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
    endfor
elseif (depth_a_len > 1);
    depth_a_values = Sim.depth_a;

    for i = 1:depth_a_len;
        close all
        Sim.horn_number = i;
        Sim.depth_a     = depth_a_values(i);
        sweep_type      = sprintf('depth_a_sweep_%.2f', depth_a_values(i));
        simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
    endfor
elseif (depth_b_len > 1);
    depth_b_values = Sim.depth_b;

    for i = 1:depth_b_len;
        close all
        Sim.horn_number = i;
        Sim.depth_b     = depth_b_values(i);
        sweep_type      = sprintf('depth_b_sweep_%.2f', depth_b_values(i));
        simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
    endfor
elseif (wg_length_len > 1);
    wg_length_values = Sim.wg_length;

    for i = 1:wg_length_len;
        close all
        Sim.horn_number = i;
        Sim.wg_length   = wg_length_values(i);
        sweep_type      = sprintf('wg_length_sweep_%.2f', wg_length_values(i));
        simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
    endfor
elseif (num_of_corrugations_len > 1);
    num_of_corrugations_values = Sim.num_of_corrugations;

    for i = 1:num_of_corrugations_len;
        close all
        Sim.horn_number         = i;
        Sim.num_of_corrugations = num_of_corrugations_values(i);
        sweep_type              = sprintf('num_of_corrugations_sweep_%d', num_of_corrugations_values(i));
        simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
    endfor
else
    close all
    Sim.horn_number = 0;
    sweep_type = 'no_sweep';
    simulate(Sim_Path, Sim_CSX, Sim, sweep_type, RUN_SIMULATION);
endif

disp(">>-------------------- Simulation fineshed! --------------------<<");
