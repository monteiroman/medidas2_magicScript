%
%
%
function [port, nf2ff] = make_horn(Sim_Path, Sim_CSX, Sim)
    %%
    %%  Documentation is still missing
    %%

    disp('>>________ Horn and Simulation Values ________<<');
    fmin    = Sim.fmin
    fmax    = Sim.fmax
    fcalc   = Sim.fcalc

    ai = Sim.ai
    bi = Sim.bi

    ao = Sim.ao
    bo = Sim.bo

    horn_number         = Sim.horn_number;

    pitch               = Sim.pitch         
    delta               = Sim.delta
    wg_length           = Sim.wg_length
    num_of_corrugations = Sim.num_of_corrugations
    straight_width      = Sim.straight_width
    cap_width           = Sim.cap_width

    exc_mode = Sim.exc_mode

    SHOW_STRUCTURE_FIGURES      = Sim.SHOW_STRUCTURE_FIGURES
    USE_CORRUGATIONS            = Sim.USE_CORRUGATIONS
    SUBSTRACT_LEFTOVERS         = Sim.SUBSTRACT_LEFTOVERS
    CHANGE_CORRUGATIONS_DEPTH   = Sim.CHANGE_CORRUGATIONS_DEPTH


    TIME_STEPS  = Sim.TIME_STEPS
    n_cell      = Sim.n_cell

    USE_PROFILE = Sim.USE_PROFILE        % 1=Linear, 2=Tangential, 3=Exponential
    output_path = Sim.output_path
    disp('>>____________________________________________<<');

    % Save variables for analysis
    config_variables_file = strcat(Sim.output_path, sprintf('/%d_horn_variables.txt',horn_number));
    save("-text", config_variables_file, "Sim");

    unit = 1e-3;                        % Units in mm
    lambda_c = 300/fcalc;               % Center frequency wavelength
    length = num_of_corrugations*pitch; % Length of horn profile
    N = length/pitch;                   % Total number of corrugations
    z = 0:pitch:length;                 % z index distance array from 0 to length of horn
    air_guard = 2;


%%_____________________________ START OF 2D FIGURES DESIGN _____________________________
    %%% Profile for A faces %%%
    switch (USE_PROFILE)
        case 1
            % Linear profile
            a_profile = ai+(ao/2-ai/2)*z/length;
        case 2
            % Tangential profile
            A = 1;
            rho = 2;
            a_profile = ai+(ao-ai)*((1-A)*(z/length)+A*power(tan((pi*z)/(4*length)),rho));
        case 3
            % Exponential profile
            a_profile = ai*exp(log(ao/ai)*(z/length));
    endswitch

    if (SHOW_STRUCTURE_FIGURES);
        structure_figure = figure('position',[600,100,900,900]);
        subplot (3, 2, 1)
        plot(z, a_profile);    
        set(gca, "linewidth",2, "fontsize", 14 )
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Linear Horn A Profile', 'FontSize', 16 );
        axis equal;
    endif

    if (USE_CORRUGATIONS);
        % Corrugations depths.
        if (CHANGE_CORRUGATIONS_DEPTH);
            d = 1:N+1;
            depth_step = (((300/fcalc)/2) - ((300/fcalc)/4)) / N;

            for i = 1:N+1;
                d(i) = ((300/fcalc)/2) - i*depth_step;
            endfor
        else
            d = zeros(1,N+1);
            d = d + (300/fcalc)/2;
        endif
    else
        d = zeros(1,N+1);
    endif

    % Generate z,y coordinates as z_for_a_profile and y_for_a_profile vector
    n = 0;
    z_for_a_profile(1) = 0;
    z_for_a_profile(2) = 0;
    for i = 1:N;
        y_for_a_profile(i+n) = a_profile(i);
        y_for_a_profile(i+n+1) = a_profile(i)+d(i);
        y_for_a_profile(i+n+2) = a_profile(i)+d(i);
        y_for_a_profile(i+n+3) = a_profile(i+1);
        y_for_a_profile(i+n+4) = a_profile(i+1);
        z_for_a_profile(i+n+2) = z_for_a_profile(i+n)+delta*pitch;
        z_for_a_profile(i+n+3) = z_for_a_profile(i+n+2);
        z_for_a_profile(i+n+4) = z_for_a_profile(i+n+3)+(1-delta)*pitch;
        z_for_a_profile(i+n+5) = z_for_a_profile(i+n+4);
        n = n+3;
    endfor

    z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
    z_for_a_profile = z_for_a_profile(1:z_number); % Truncate z axis data points to equal y_for_a_profile vector length

    if (SHOW_STRUCTURE_FIGURES);
        subplot (3, 2, 3)
        plot(z_for_a_profile,y_for_a_profile);
        set(gca, "linewidth",2, "fontsize", 14 )
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Corrugation A Profile', 'FontSize', 16 );
        axis equal;   % Scale axis equally for aspect ratio 1:1
    endif

    % Add the rest of the geometry to create a closed path
    % a_offset is the inner horn profile shifted up to give the horn a thickness
    a_offset = a_profile.+(lambda_c/2+2);
    % figure;                    % Uncomment these three lines for debugging
    % plot(z, a_offset);
    % axis equal;

    a_volume_y_end = a_offset(end);

    % Add vertical surface at horn aperture
    z_for_a_profile = [z_for_a_profile, z_for_a_profile(z_number)];
    y_for_a_profile = [y_for_a_profile, a_offset(N)];
    mesh_a = y_for_a_profile;                 % radmesh to fix mesh lines to corrugations
    % figure;                    % Uncomment these three lines for debugging
    % plot(z_for_a_profile, y_for_a_profile);
    % axis equal;

    % Flip outer surface profile so that widest horn dimensions comes next in the outline coordinates
    outer_surface = fliplr(a_offset);
    z_flip = fliplr(z);
    extent = z_for_a_profile(end);  % Fudge to make horn aperture planar for ring loaded slot MC
    z_flip(1) = extent; % Fudge to make horn aperture planar for ring loaded slot MC
    % Add outer profile and waveguide to horn
    z_for_a_profile = [z_for_a_profile, z_flip,         -wg_length,               -wg_length,   0];
    y_for_a_profile = [y_for_a_profile, outer_surface,  ai+(lambda_c/2+2),  ai,     ai];

    if (SHOW_STRUCTURE_FIGURES);
        subplot (3, 2, 5)
        plot(z_for_a_profile,y_for_a_profile);
        set(gca, "linewidth",2, "fontsize", 14 )
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Complete Corrugated Horn A Profile', 'FontSize', 16 );
        axis equal;   % Scale axis equally for aspect ratio 1:1
    endif

    % Substraction volume for A
    z_for_a_subs_profile = [z_flip,         -wg_length-air_guard,   -wg_length-air_guard,       0];  
    z_for_a_subs_profile = [z_for_a_subs_profile,   z_flip(1),                  z_flip(1)];
    y_for_a_subs_profile = [outer_surface,  ai+(lambda_c/2+2),      outer_surface(1)+air_guard, outer_surface(1)+air_guard]; 
    y_for_a_subs_profile = [y_for_a_subs_profile,   outer_surface(1)+air_guard, outer_surface(1)];

    %%% Profile for B faces %%%
    switch (USE_PROFILE)
        case 1
            % Linear profile
            b_profile = bi+(bo/2-bi/2)*z/length;
        case 2
            % Tangential profile
            A = 1;
            rho = 2;
            b_profile = bi+(bo-bi)*((1-A)*(z/length)+A*power(tan((pi*z)/(4*length)),rho));
        case 3
            % Exponential profile
            b_profile = bi*exp(log(bo/bi)*(z/length));
    endswitch

    if (SHOW_STRUCTURE_FIGURES);
        subplot (3, 2, 2)
        plot(z, b_profile);    
        set(gca, "linewidth",2, "fontsize", 14 )
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Linear Horn B Profile', 'FontSize', 16 );
        axis equal;
    endif

    if (USE_CORRUGATIONS);
        % Corrugations depths.
        if (CHANGE_CORRUGATIONS_DEPTH);
            d = 1:N+1;
            depth_step = (((300/fcalc)/2) - ((300/fcalc)/4)) / N;

            for i = 1:N+1;
                d(i) = ((300/fcalc)/2) - i*depth_step;
            endfor
        else
            d = zeros(1,N+1);
            d = d + (300/fcalc)/2;
        endif
    else
        d = zeros(1,N+1);
    endif

    % Generate z,y coordinates as z_for_b_profile and y_for_b_profile vector
    n = 0;
    z_for_b_profile(1) = 0;
    z_for_b_profile(2) = 0;
    for i = 1:N;
        y_for_b_profile(i+n) = b_profile(i);
        y_for_b_profile(i+n+1) = b_profile(i)+d(i);
        y_for_b_profile(i+n+2) = b_profile(i)+d(i);
        y_for_b_profile(i+n+3) = b_profile(i+1);
        y_for_b_profile(i+n+4) = b_profile(i+1);
        z_for_b_profile(i+n+2) = z_for_b_profile(i+n)+delta*pitch;
        z_for_b_profile(i+n+3) = z_for_b_profile(i+n+2);
        z_for_b_profile(i+n+4) = z_for_b_profile(i+n+3)+(1-delta)*pitch;
        z_for_b_profile(i+n+5) = z_for_b_profile(i+n+4);
        n = n+3;
    endfor

    z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
    z_for_b_profile = z_for_b_profile(1:z_number); % Truncate z axis data points to equal y_for_b_profile vector length

    if (SHOW_STRUCTURE_FIGURES);
        subplot (3, 2, 4)
        plot(z_for_b_profile,y_for_b_profile);
        set(gca, "linewidth",2, "fontsize", 14 )
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Corrugation B Profile', 'FontSize', 16 );
        axis equal;   % Scale axis equally for aspect ratio 1:1
    endif

    % Add the rest of the geometry to create a closed path
    % b_offset is the inner horn profile shifted up to give the horn a thickness
    b_offset = b_profile.+(lambda_c/2+2);
    % figure;                    % Uncomment these three lines for debugging
    % plot(z, b_offset);
    % axis equal;

    b_volume_y_end = b_offset(end);

    % Add vertical surface at horn aperture
    z_for_b_profile = [z_for_b_profile, z_for_b_profile(z_number)];
    y_for_b_profile = [y_for_b_profile, b_offset(N)];
    mesh_b = y_for_b_profile;                 % radmesh to fix mesh lines to corrugations
    % figure;                    % Uncomment these three lines for debugging
    % plot(z_for_b_profile, y_for_b_profile);
    % axis equal;

    % Flip outer surface profile so that widest horn dimensions comes next in the outline coordinates
    outer_surface = fliplr(b_offset);
    z_flip = fliplr(z);
    extent = z_for_b_profile(end);  % Fudge to make horn aperture planar for ring loaded slot MC
    z_flip(1) = extent; % Fudge to make horn aperture planar for ring loaded slot MC
    % Add outer profile and waveguide to horn
    z_for_b_profile = [z_for_b_profile, z_flip,         -wg_length,               -wg_length,   0];
    y_for_b_profile = [y_for_b_profile, outer_surface,  bi+(lambda_c/2+2),  bi,     bi];

    if (SHOW_STRUCTURE_FIGURES);
        subplot (3, 2, 6)
        plot(z_for_b_profile,y_for_b_profile);
        set(gca, "linewidth",2, "fontsize", 14 )
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Complete Corrugated Horn B Profile', 'FontSize', 16 );
        axis equal;   % Scale axis equally for aspect ratio 1:1

        structure_output = strcat(Sim.output_path, sprintf('/%d_structure_output.svg',horn_number));
        print (structure_figure, structure_output, "-dsvg", "-Sxsize=900");
    endif

    % Substraction volume for B
    z_for_b_subs_profile = [z_flip,         -wg_length-air_guard,   -wg_length-air_guard,       0];  
    y_for_b_subs_profile = [outer_surface,  bi+(lambda_c/2+2),      outer_surface(1)+air_guard, outer_surface(1)+air_guard];
    z_for_b_subs_profile = [z_for_b_subs_profile,   z_flip(1),                  z_flip(1)];
    y_for_b_subs_profile = [y_for_b_subs_profile,   outer_surface(1)+air_guard, outer_surface(1)];

    %% Generate end cap to prevent the radiation coming out of the back of the horn
    % Extract the end point of the structure
    cap_point_z_1 = [z_flip, -wg_length](end);
    cap_point_y_1 = [outer_surface, ai+(lambda_c/2+2)](end) - ai/2;
    % From here we have to move in "y" to to met the "y" end of the cap.
    cap_point_z_2 = cap_point_z_1;
    cap_point_y_2 = -cap_point_y_1;
    % We define here the width of the cap
    cap_point_z_3 = cap_point_z_1 + cap_width;
    cap_point_y_3 = cap_point_y_1;
    cap_point_z_4 = cap_point_z_2 + cap_width;
    cap_point_y_4 = cap_point_y_2;

    cap_y = [cap_point_y_1, cap_point_y_2, cap_point_y_4, cap_point_y_3, cap_point_y_1];
    cap_z = [cap_point_z_1, cap_point_z_2, cap_point_z_4, cap_point_z_3, cap_point_z_1];

    %%%_____________________________ END OF 2D FIGURES DESIGN _____________________________


    % openEMS setup begins here
    % EM related physical constants
    physical_constants;

    % frequency range of interest
    f_start =  fmin*1e9;
    f_stop  =  fmax*1e9;

    % frequency to calculate fields
    f0 = fcalc*1e9;

    %% setup FDTD parameter & excitation function
    FDTD = InitFDTD( 'NrTS', TIME_STEPS, 'EndCriteria', 0.5e-3 );
    FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
    BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % FDTD Boundary Conditions:
                                                                    % http://openems.de/index.php/FDTD_Boundary_Conditions
    FDTD = SetBoundaryCond(FDTD, BC);

    %% setup CSXCAD geometry & mesh
    max_res = c0/(f_stop)/unit/n_cell;  % cell size: lambda/20
    CSX = InitCSX();                    % Initialise CSX structure

    % Calculate lambda/4 at lowest frequency to use as distance to nf2ff surfaces
    lambda_max = c0/f_start/unit/4;

    % Create fixed lines for the simulation box, structure and port
    % Info at this link: http://openems.de/index.php/FDTD_Mesh
    mesh.y = [(-b_offset(end)-(9*max_res)-lambda_max) -mesh_b(1:4:end)+bi/2 0 mesh_b(1:4:end)-bi/2 (b_offset(end)+(9*max_res)+
                                                                                                            lambda_max)];
    mesh.y = SmoothMeshLines( mesh.y, max_res, 1.5); % Create a smooth mesh between specified fixed mesh lines

    mesh.x = [(-a_offset(end)-(9*max_res)-lambda_max) -mesh_a(1:4:end)+ai/2 0 mesh_a(1:4:end)-ai/2 (a_offset(end)+(9*max_res)+
                                                                                                            lambda_max)];
    mesh.x = SmoothMeshLines( mesh.x, max_res, 1.5); % Create a smooth mesh between specified fixed mesh lines

    % % Create fixed lines for the simulation box,port and given number of lines inside the horn
    mesh.z = [-wg_length-lambda_max-(9*max_res) -wg_length-1 -wg_length -wg_length+10 0 z_for_a_profile(1:2:z_number) length+2*lambda_max+(9*max_res)];
    mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );



    % mesh.y = [(-b_offset(end)-(9*max_res)-lambda_max) -bi/2 0 bi/2 (b_offset(end)+(9*max_res)+lambda_max)];
    % mesh.y = SmoothMeshLines( mesh.y, max_res, 1.5); % Create a smooth mesh between specified fixed mesh lines

    % mesh.x = [(-a_offset(end)-(9*max_res)-lambda_max) -ai/2 0 ai/2 (a_offset(end)+(9*max_res)+lambda_max)];
    % mesh.x = SmoothMeshLines( mesh.x, max_res, 1.5); % Create a smooth mesh between specified fixed mesh lines

    % % Create fixed lines for the simulation box,port and given number of lines inside the horn
    % mesh.z = [-wg_length-lambda_max-(9*max_res) 0 length+2*lambda_max+(9*max_res)];
    % mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );






    CSX = DefineRectGrid( CSX, unit, mesh );

    %% ----->> Create Horn Geometry <<-----
    %% Horn + waveguide
    CSX = AddMetal(CSX, 'Corrugated_Horn');
    CSX = AddMetal(CSX, 'Cap');
    if (SUBSTRACT_LEFTOVERS);
        % openEMSs way to subtract volumes: https://openems.de/index.php/Metal_sheet_with_cylindrical_holes.html
        %                       primitives: https://openems.de/index.php/Primitives.html#Coordinate_System_Definition
        CSX = AddMaterial( CSX, 'Air' );
        CSX = SetMaterialProperty( CSX, 'Air', 'Epsilon', 1, 'Mue', 1 );
    endif

    corrugated_coords_a = [y_for_a_profile; z_for_a_profile];
    substract_coords_a  = [y_for_a_subs_profile; z_for_a_subs_profile];
    corrugated_coords_b = [y_for_b_profile; z_for_b_profile];
    substract_coords_b  = [y_for_b_subs_profile; z_for_b_subs_profile];
    guard = 2;

    % Corrugated A walls
    % https://openems.de/index.php/Polygon.html
    CSX = AddLinPoly( CSX, 'Corrugated_Horn', 5, 1, 0, corrugated_coords_a, b_volume_y_end*2-bi, 'Transform', {
        'Rotate_Y', -pi/2, 'Translate',[ num2str(ai/2) ',' num2str(-(b_volume_y_end*2-bi)/2) ',0']
        });
    CSX = AddLinPoly( CSX, 'Corrugated_Horn', 5, 1, 0, corrugated_coords_a, b_volume_y_end*2-bi, 'Transform', {
        'Rotate_Y', pi/2, 'Rotate_X', pi, 'Translate',[ num2str(-ai/2) ',' num2str((b_volume_y_end*2-bi)/2) ',0']
        });

    % Subtract horn A walls outside geometry
    if (SUBSTRACT_LEFTOVERS);
        CSX = AddLinPoly( CSX, 'Air', 10, 1, 0, substract_coords_a, b_volume_y_end*2-bi, 'Transform', {
            'Rotate_Y', -pi/2, 'Translate',[ num2str(ai/2) ',' num2str(-(b_volume_y_end*2-bi)/2) ',0']
            });
        CSX = AddLinPoly( CSX, 'Air', 10, 1, 0, substract_coords_a, b_volume_y_end*2-bi, 'Transform', {
            'Rotate_Y', pi/2, 'Rotate_X', pi, 'Translate',[ num2str(-ai/2) ',' num2str((b_volume_y_end*2-bi)/2) ',0']
            });
    endif

    % Corrugated B walls
    CSX = AddLinPoly( CSX, 'Corrugated_Horn', 5, 1, 0, corrugated_coords_b, a_volume_y_end*2-ai, 'Transform', {
        'Rotate_Y', -pi/2, 'Rotate_Z', -pi/2, 'Translate',[ num2str(-(a_volume_y_end*2-ai)/2) ',' num2str(-bi/2) ',0']
        });
    CSX = AddLinPoly( CSX, 'Corrugated_Horn', 5, 1, 0, corrugated_coords_b, a_volume_y_end*2-ai, 'Transform', {
        'Rotate_Y', -pi/2, 'Rotate_Z', pi/2, 'Translate',[ num2str((a_volume_y_end*2-ai)/2) ',' num2str(bi/2) ',0']
        });

    % Subtract horn B walls outside geometry
    if (SUBSTRACT_LEFTOVERS);
        CSX = AddLinPoly( CSX, 'Air', 10, 1, 0, substract_coords_b, a_volume_y_end*2-ai, 'Transform', {
            'Rotate_Y', -pi/2, 'Rotate_Z', -pi/2, 'Translate',[ num2str(-(a_volume_y_end*2-ai)/2) ',' num2str(-bi/2) ',0']
            });
        CSX = AddLinPoly( CSX, 'Air', 10, 1, 0, substract_coords_b, a_volume_y_end*2-ai, 'Transform', {
            'Rotate_Y', -pi/2, 'Rotate_Z', pi/2, 'Translate',[ num2str((a_volume_y_end*2-ai)/2) ',' num2str(bi/2) ',0']
            });
    endif


    %% End cap to prevent the radiation coming out of the back of the horn
    CSX = AddMetal(CSX, 'Cap');
    cap_coords = [cap_y; cap_z];

    CSX = AddLinPoly( CSX, 'Cap', 10, 1, 0, cap_coords, bi + 2 * straight_width, 'Transform',{
        'Rotate_Y', -pi/2, 'Translate',[ '0,' num2str(-bi/2 - straight_width) ',' num2str(-cap_width)]
        });
    %% ----->> End of model geometry <<-----


    %% ----->> Excitation, dumpboxes and nf2ff <<-----
    % Apply the excitation
    start=[-ai/2 -bi/2 -wg_length];
    stop =[ai/2 bi/2 -wg_length+10];
    %% Ports information at link:
    %       https://openems.de/index.php/Ports.html#Rectangular_Waveguide_Ports
    [CSX, port] = AddRectWaveGuidePort( CSX, 0, 1, start, stop, 'z', ai*unit, bo*unit, exc_mode, 1);

    % Dump box for Electric field at Phi=0 (vertical cut)
    CSX = AddDump(CSX,'Et_V_dump', 'SubSampling', '4,4,4');
    start=[0 (-a_offset(end)-lambda_max) (-wg_length-lambda_max)];
    stop =[0 (a_offset(end)+lambda_max) (length+2*lambda_max)];
    CSX = AddBox(CSX,'Et_V_dump',0,start,stop);

    % Dump box for Electric field at Phi=90 (horizontal cut)
    CSX = AddDump(CSX,'Et_H_dump', 'SubSampling', '4,4,4');
    start=[(-a_offset(end)-lambda_max) 0 (-wg_length-lambda_max)];
    stop =[(a_offset(end)+lambda_max) 0 (length+2*lambda_max)];
    CSX = AddBox(CSX,'Et_H_dump',0,start,stop);

    % nf2ff calc
    start = [mesh.x(9) mesh.y(9) mesh.z(9)];
    stop  = [mesh.x(end-8) mesh.y(end-8) mesh.z(end-8)];
    [CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'Directions', [1 1 1 1 1 1], 'OptResolution', max_res*4);

    %% ----->> End of excitation, dumpboxes and nf2ff definitions <<-----

    %% ----->> Prepare simulation folder <<----- 
    confirm_recursive_rmdir(0);                                 % No ask if removedirectory
    [status, message, messageid] = rmdir( Sim_Path, 's');       % Clear previous directory
    [status, message, messageid] = mkdir( Sim_Path );           % Create empty simulation folder

    % Write openEMS compatible xml-file
    WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

    % Show structure
    CSXGeomPlot([Sim_Path '/' Sim_CSX], ['--export-STL=tmp']);

    %% ----->> End of simulation folder <<----- 
endfunction
