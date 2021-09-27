%% Clean workspace
close all
clear
clc

addpath('src');

%%
%%
%%
function simulate(Sim_Path, Sim_CSX, Sim, RUN_SIMULATION)
    %%
    %%  Documentation is still missing
    %%

    %% ----->> Output simulation particular folder <<----- 
    Sim.output_path = strcat(Sim.output_path, sprintf('/Corrugated_Horn_%d',Sim.horn_number));
    [status, message, messageid] = mkdir(Sim.output_path);

    %% ----->> Make horn structure and save openEMS XML simulation file at "<Sim_Path>/<Sim_CSX>"  <<-----
    [port, nf2ff] = make_horn(Sim_Path, Sim_CSX, Sim);

    %% ----->> Simulate and store 3D files <<-----
    if (RUN_SIMULATION);
        run_simulation(Sim_Path, Sim_CSX, Sim, port, nf2ff);
        movefile(strcat(Sim_Path, '/*.stl'), Sim.output_path);
        movefile(strcat(Sim_Path, '/*.vtk'), Sim.output_path);
    endif
endfunction

%------------------------------------------------------------------------------%
%                                   MagicScript                                %
%------------------------------------------------------------------------------%
%%
%%  Documentation is still missing
%%
RUN_SIMULATION = 1;
PLOT_OUTPUT_SAME_WINDOW = 0;


output_path = 'outputs';
Sim_Path    = 'tmp';
Sim_CSX     = 'Corrugated_Horn.xml';

Sim.fmin    = 8;                     % Minimum frequency in GHz
Sim.fmax    = 12;                    % Maximum frequency in GHz
Sim.fcalc   = 10;                    % Frequency to calculate fields

Sim.ai      = 22.86;
Sim.bi      = 10.16;

Sim.ao      = [125:10:225];
Sim.bo      = 120;

Sim.horn_number = 1;

Sim.pitch               = 4;         
Sim.delta               = 0.75;      % Pitch to width ratio
Sim.wg_length           = 60;       % Length of feeding waveguide
Sim.num_of_corrugations = 40;
Sim.straight_width      = 2;
Sim.cap_width           = 2;

Sim.exc_mode = 'TE10';

Sim.SHOW_STRUCTURE_FIGURES      = 1;
Sim.USE_CORRUGATIONS            = 1;
Sim.SUBSTRACT_LEFTOVERS         = 1;
Sim.CHANGE_CORRUGATIONS_DEPTH   = 1;

Sim.TIME_STEPS  = 5000;
Sim.n_cell      = 20;                    % cell size: lambda/n_cell

Sim.USE_PROFILE = 1;                     % 1=Linear, 2=Tangential, 3=Exponential
Sim.output_path = output_path;


%% ----->> Output simulation generic folder <<----- 
confirm_recursive_rmdir(0);                                 % No ask if removedirectory
[status, message, messageid] = rmdir( output_path, 's');    % Clear previous directory
[status, message, messageid] = mkdir( output_path );        % Create empty simulation folder

ai_len = length(Sim.ai);
bi_len = length(Sim.bi);
ao_len = length(Sim.ao);
bo_len = length(Sim.bo);

if (ai_len > 1);
    ai_values = Sim.ai;

    for i = 1:ai_len+1;
        Sim.horn_number = i;
        Sim.ai = ai_values(i);
        simulate(Sim_Path, Sim_CSX, Sim, RUN_SIMULATION);
    endfor
elseif (bi_len > 1);
    bi_values = Sim.bi;

    for i = 1:bi_len+1;
        Sim.horn_number = i;
        Sim.bi = bi_values(i);
        simulate(Sim_Path, Sim_CSX, Sim, RUN_SIMULATION);
    endfor
elseif (ao_len > 1);
    ao_values = Sim.ao;

    for i = 1:ao_len+1;
        Sim.horn_number = i;
        Sim.ao = ao_values(i);
        simulate(Sim_Path, Sim_CSX, Sim, RUN_SIMULATION);
    endfor
elseif (bo_len > 1);
    bo_values = Sim.bo;

    for i = 1:bo_len+1;
        Sim.horn_number = i;
        Sim.bo = bo_values(i);
        simulate(Sim_Path, Sim_CSX, Sim, RUN_SIMULATION);
    endfor
endif

disp(">>-------------------- Sweep fineshed! --------------------<<");
