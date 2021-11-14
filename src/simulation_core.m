%%
%%  Author: Tiago Monteiro
%%
%% 
function simulation_core(Sim, RUN_SIMULATION)
    %%  
    %%  Simulation sweep function. Checks parameters of interest lengths and loops its values if
    %%  there are more than one.
    %%  
    %%  Parameters:
    %%          Sim: Simulation structures parameters
    %%          RUN_SIMULATION: Flag that determines if the simulation must be 
    %%                          done.  
    %%  
    %%  Returns:
    %%          void.
    %%

    %% ----->> Generic simulation output folder <<----- 
    confirm_recursive_rmdir(0);                                     % Do not asks if remove directory
    [status, message, messageid] = rmdir( Sim.output_path, 's');    % Clear previous directory
    [status, message, messageid] = mkdir( Sim.output_path );        % Create empty simulation folder

    if (Sim.MAKE_HORN);
        %% ----->> Check parameters and sweep if necessary <<----- 
        ai_len = length(Sim.horn_ai);
        bi_len = length(Sim.horn_bi);
        ao_len = length(Sim.horn_ao);
        bo_len = length(Sim.horn_bo);
        corr_step_len       = length(Sim.horn_corr_step);
        delta_len           = length(Sim.horn_delta);
        depth_a_len         = length(Sim.horn_depth_a);
        a_jump_len          = length(Sim.horn_a_jump);
        A_wall_tan_A_len    = length(Sim.horn_A_wall_tan_A);
        A_wall_tan_rho_len  = length(Sim.horn_A_wall_tan_rho);
        first_a_length_len  = length(Sim.horn_first_a_length);
        first_ao_len        = length(Sim.horn_first_ao);
        depth_b_len         = length(Sim.horn_depth_b);
        b_jump_len          = length(Sim.horn_b_jump);
        B_wall_tan_A_len    = length(Sim.horn_B_wall_tan_A);
        B_wall_tan_rho_len  = length(Sim.horn_B_wall_tan_rho);
        first_b_length_len  = length(Sim.horn_first_b_length);
        first_bo_len        = length(Sim.horn_first_bo);
        wg_length_len       = length(Sim.horn_wg_length);
        num_of_corrugations_len = length(Sim.horn_num_of_corrugations);

        if (ai_len > 1);
            ai_values = Sim.horn_ai;

            for i = 1:ai_len;
                close all                                   % Prevent memory leakage
                Sim.horn_number = i;
                Sim.horn_ai     = ai_values(i);
                sweep_type      = sprintf('ai_sweep_%.2f', ai_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (bi_len > 1);
            bi_values = Sim.horn_bi;

            for i = 1:bi_len;
                close all
                Sim.horn_number = i;
                Sim.horn_bi     = bi_values(i);
                sweep_type      = sprintf('bi_sweep_%.2f', bi_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (ao_len > 1);
            ao_values = Sim.horn_ao;

            for i = 1:ao_len;
                close all
                Sim.horn_number = i;
                Sim.horn_ao     = ao_values(i);
                sweep_type      = sprintf('ao_sweep_%.2f', ao_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (bo_len > 1);
            bo_values = Sim.horn_bo;

            for i = 1:bo_len;
                close all
                Sim.horn_number = i;
                Sim.horn_bo     = bo_values(i);
                sweep_type      = sprintf('bo_sweep_%.2f', bo_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (corr_step_len > 1);
            corr_step_values = Sim.horn_corr_step;

            for i = 1:corr_step_len;
                close all
                Sim.horn_number     = i;
                Sim.horn_corr_step  = corr_step_values(i);
                sweep_type          = sprintf('corr_step_sweep_%.2f', corr_step_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (delta_len > 1);
            delta_values = Sim.horn_delta;

            for i = 1:delta_len;
                close all
                Sim.horn_number = i;
                Sim.horn_delta  = delta_values(i);
                sweep_type      = sprintf('delta_sweep_%.2f', delta_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (depth_a_len > 1);
            depth_a_values = Sim.horn_depth_a;

            for i = 1:depth_a_len;
                close all
                Sim.horn_number  = i;
                Sim.horn_depth_a = depth_a_values(i);
                sweep_type       = sprintf('depth_a_sweep_%.2f', depth_a_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (a_jump_len > 1);
            a_jump_values = Sim.horn_a_jump;

            for i = 1:a_jump_len;
                close all
                Sim.horn_number = i;
                Sim.horn_a_jump = a_jump_values(i);
                sweep_type      = sprintf('a_jump_sweep_%.2f', a_jump_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (A_wall_tan_A_len > 1);
            A_wall_tan_A_values = Sim.horn_A_wall_tan_A;

            for i = 1:A_wall_tan_A_len;
                close all
                Sim.horn_number         = i;
                Sim.horn_A_wall_tan_A   = A_wall_tan_A_values(i);
                sweep_type              = sprintf('A_wall_tan_A_sweep_%.2f', A_wall_tan_A_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (A_wall_tan_rho_len > 1);
            A_wall_tan_rho_values = Sim.horn_A_wall_tan_rho;

            for i = 1:A_wall_tan_rho_len;
                close all
                Sim.horn_number         = i;
                Sim.horn_A_wall_tan_rho = A_wall_tan_rho_values(i);
                sweep_type              = sprintf('A_wall_tan_rho_sweep_%.2f', A_wall_tan_rho_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (first_a_length_len > 1);
            first_a_length_values = Sim.horn_first_a_length;

            for i = 1:first_a_length_len;
                close all
                Sim.horn_number         = i;
                Sim.horn_first_a_length = first_a_length_values(i);
                sweep_type              = sprintf('first_a_length_sweep_%.2f', first_a_length_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (first_ao_len > 1);
            first_ao_values = Sim.horn_first_ao;

            for i = 1:first_ao_len;
                close all
                Sim.horn_number     = i;
                Sim.horn_first_ao   = first_ao_values(i);
                sweep_type          = sprintf('first_ao_sweep_%.2f', first_ao_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (depth_b_len > 1);
            depth_b_values = Sim.horn_depth_b;

            for i = 1:depth_b_len;
                close all
                Sim.horn_number     = i;
                Sim.horn_depth_b    = depth_b_values(i);
                sweep_type          = sprintf('depth_b_sweep_%.2f', depth_b_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (b_jump_len > 1);
            b_jump_values = Sim.horn_b_jump;

            for i = 1:b_jump_len;
                close all
                Sim.horn_number = i;
                Sim.horn_b_jump = b_jump_values(i);
                sweep_type      = sprintf('b_jump_sweep_%.2f', b_jump_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (B_wall_tan_A_len > 1);
            B_wall_tan_A_values = Sim.horn_B_wall_tan_A;

            for i = 1:B_wall_tan_A_len;
                close all
                Sim.horn_number         = i;
                Sim.horn_B_wall_tan_A   = B_wall_tan_A_values(i);
                sweep_type              = sprintf('B_wall_tan_A_sweep_%.2f', B_wall_tan_A_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (B_wall_tan_rho_len > 1);
            B_wall_tan_rho_values = Sim.horn_B_wall_tan_rho;

            for i = 1:B_wall_tan_rho_len;
                close all
                Sim.horn_number         = i;
                Sim.horn_B_wall_tan_rho = B_wall_tan_rho_values(i);
                sweep_type              = sprintf('B_wall_tan_rho_sweep_%.2f', B_wall_tan_rho_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (first_b_length_len > 1);
            first_b_length_values = Sim.horn_first_b_length;

            for i = 1:first_b_length_len;
                close all
                Sim.horn_number         = i;
                Sim.horn_first_b_length = first_b_length_values(i);
                sweep_type              = sprintf('first_b_length_sweep_%.2f', first_b_length_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (first_bo_len > 1);
            first_bo_values = Sim.horn_first_bo;

            for i = 1:first_bo_len;
                close all
                Sim.horn_number     = i;
                Sim.horn_first_bo   = first_bo_values(i);
                sweep_type          = sprintf('first_bo_sweep_%.2f', first_bo_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (wg_length_len > 1);
            wg_length_values = Sim.horn_wg_length;

            for i = 1:wg_length_len;
                close all
                Sim.horn_number     = i;
                Sim.horn_wg_length  = wg_length_values(i);
                sweep_type          = sprintf('wg_length_sweep_%.2f', wg_length_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (num_of_corrugations_len > 1);
            num_of_corrugations_values = Sim.horn_num_of_corrugations;

            for i = 1:num_of_corrugations_len;
                close all
                Sim.horn_number              = i;
                Sim.horn_num_of_corrugations = num_of_corrugations_values(i);
                sweep_type                   = sprintf('num_of_corrugations_sweep_%d', num_of_corrugations_values(i));
                simulate_horn(Sim, sweep_type, RUN_SIMULATION);
            endfor
        else
            close all
            Sim.horn_number = 0;
            sweep_type = 'no_sweep';
            simulate_horn(Sim, sweep_type, RUN_SIMULATION);
        endif
    endif

    if (RUN_SIMULATION == 0)  % Gives AppCSXCad time to open the horn model.
        pause(5);
    endif

    if (Sim.MAKE_ADAPT);
        %% ----->> Check parameters and sweep if necessary <<----- 
        BackShort_len = length(Sim.adapt_BackShort);
        ProbeDepth_len = length(Sim.adapt_ProbeDepth);
        adapt_a_len = length(Sim.adapt_a);
        adapt_b_len = length(Sim.adapt_b);
        adapt_dielectric_intrusion_len = length(Sim.adapt_dielectric_intrusion);
        adapt_ProbeRad_len = length(Sim.adapt_ProbeRad);

        if (BackShort_len > 1);
            BackShort_values = Sim.adapt_BackShort;

            for i = 1:BackShort_len;
                close all                                   % Prevent memory leakage
                Sim.adapt_number    = i;
                Sim.adapt_BackShort = BackShort_values(i);
                sweep_type          = sprintf('BackShort_sweep_%.2f', BackShort_values(i));
                simulate_adapter(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (ProbeDepth_len >1);
            ProbeDepth_values = Sim.adapt_ProbeDepth;

            for i = 1:ProbeDepth_len;
                close all                                   
                Sim.adapt_number     = i;
                Sim.adapt_ProbeDepth = ProbeDepth_values(i);
                sweep_type           = sprintf('ProbeDepth_sweep_%.2f', ProbeDepth_values(i));
                simulate_adapter(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (adapt_a_len >1);
            adapt_a_values = Sim.adapt_a;

            for i = 1:adapt_a_len;
                close all                                   
                Sim.adapt_number = i;
                Sim.adapt_a      = adapt_a_values(i);
                sweep_type       = sprintf('adapt_a_sweep_%.2f', adapt_a_values(i));
                simulate_adapter(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (adapt_b_len >1);
            adapt_b_values = Sim.adapt_b;

            for i = 1:adapt_b_len;
                close all                                   
                Sim.adapt_number = i;
                Sim.adapt_b      = adapt_b_values(i);
                sweep_type       = sprintf('adapt_b_sweep_%.2f', adapt_b_values(i));
                simulate_adapter(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (adapt_dielectric_intrusion_len >1);
            adapt_dielectric_intrusion_values = Sim.adapt_dielectric_intrusion;

            for i = 1:adapt_dielectric_intrusion_len;
                close all                                   
                Sim.adapt_number                = i;
                Sim.adapt_dielectric_intrusion  = adapt_dielectric_intrusion_values(i);
                sweep_type                      = sprintf('adapt_dielectric_intrusion_sweep_%.2f', adapt_dielectric_intrusion_values(i));
                simulate_adapter(Sim, sweep_type, RUN_SIMULATION);
            endfor
        elseif (adapt_ProbeRad_len >1);
            adapt_ProbeRad_values = Sim.adapt_ProbeRad;

            for i = 1:adapt_ProbeRad_len;
                close all                                   
                Sim.adapt_number   = i;
                Sim.adapt_ProbeRad = adapt_ProbeRad_values(i);
                sweep_type         = sprintf('adapt_ProbeRad_sweep_%.2f', adapt_ProbeRad_values(i));
                simulate_adapter(Sim, sweep_type, RUN_SIMULATION);
            endfor
        else
            close all
            Sim.horn_number = 0;
            sweep_type = 'no_sweep';
            simulate_adapter(Sim, sweep_type, RUN_SIMULATION);
        endif
    endif
    
    pause(5);  % Gives AppCSXCad time to open the model.
endfunction

%%
%%  Simulation function
%%
function simulate_horn(Sim, sweep_type, RUN_SIMULATION)
    %%  
    %%  This function makes the simulation directory, calls make_horn for make 
    %%  the 3D structure and then calls run_horn_simulation if it is required to. 
    %%  
    %%  Parameters:
    %%          Sim: Horn simulation parameters.
    %%          sweep_type: String detailing sweep type.
    %%          RUN_SIMULATION: Flag that determines if the simulation must be 
    %%                          done.
    %%  Returns:
    %%          void
    %%

    %% ----->> Output simulation particular folder <<-----
    if (Sim.horn_number == 0);
        sweep_path = strcat(sweep_type, '_Horn_output');
    else
        sweep_path = strcat(sweep_type, sprintf('_Horn_n%d_output',Sim.horn_number));
    endif
    Sim.output_path = strcat(Sim.output_path, sweep_path);
    [status, message, messageid] = mkdir(Sim.output_path);
    [status, message, messageid] = mkdir(strcat(Sim.output_path, '/svg_plots'));
    [status, message, messageid] = mkdir(strcat(Sim.output_path, '/jpg_plots'));

    %% ----->> Make horn structure and save openEMS XML simulation file at "<Sim_Path>/<Sim_CSX>"  <<-----
    [port, nf2ff] = make_horn(Sim);

    %% ----->> Simulate and store 3D files <<-----
    if (RUN_SIMULATION);
        run_horn_simulation(Sim, port, nf2ff);

        % Move output files to output folder
        movefile(strcat(Sim.Sim_Path, '/*.txt'), Sim.output_path);
        disp("--->> Moving 3D files to output folder... <<---");
        out_3d = strcat(Sim.output_path, '/3D_Outputs');
        [status, message, messageid] = mkdir(out_3d);
        movefile(strcat(Sim.Sim_Path, '/*.stl'), out_3d);
        movefile(strcat(Sim.Sim_Path, '/*.vtk'), out_3d);
    endif
endfunction

%%
%%  Simulation function
%%
function simulate_adapter(Sim, sweep_type, RUN_SIMULATION)
    %%  
    %%  This function makes the simulation directory, calls make_adapter for make 
    %%  the 3D structure and then calls run_adapter_simulation if it is required to. 
    %%  
    %%  Parameters:
    %%          Sim: Adapter simulation parameters.
    %%          sweep_type: String detailing sweep type.
    %%          RUN_SIMULATION: Flag that determines if the simulation must be 
    %%                          done.
    %%  Returns:
    %%          void
    %%

    %% ----->> Output simulation particular folder <<-----
    if (Sim.adapt_number == 0);
        sweep_path = strcat(sweep_type, '_Adapter_output');
    else
        sweep_path = strcat(sweep_type, sprintf('_Adapter_n%d_output',Sim.adapt_number));
    endif
    Sim.output_path = strcat(Sim.output_path, sweep_path);
    [status, message, messageid] = mkdir(Sim.output_path);
    [status, message, messageid] = mkdir(strcat(Sim.output_path, '/svg_plots'));
    [status, message, messageid] = mkdir(strcat(Sim.output_path, '/jpg_plots'));

    %% ----->> Make horn structure and save openEMS XML simulation file at "<Sim_Path>/<Sim_CSX>"  <<-----
    [port, freq] = make_adapter(Sim);

    %% ----->> Simulate and store 3D files <<-----
    if (RUN_SIMULATION);
        run_adapter_simulation(Sim, port, freq);

        % Move output files to output folder
        movefile(strcat(Sim.Sim_Path, '/*.txt'), Sim.output_path);
        disp("--->> Moving 3D files to output folder... <<---");
        out_3d = strcat(Sim.output_path, '/3D_Outputs');
        [status, message, messageid] = mkdir(out_3d);
        movefile(strcat(Sim.Sim_Path, '/*.stl'), out_3d);
    endif
endfunction