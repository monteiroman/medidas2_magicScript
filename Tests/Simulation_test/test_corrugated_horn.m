close all
clear
clc
% https://www.miwv.com/x-band-horn-antennas-wr-90-8-2-12-4-ghz/


%%_____________________________ USER EDITABLE PARAMETERS _____________________________
fmin    = 8                     % Minimum frequency in GHz
fmax    = 12                    % Maximum frequency in GHz
fcalc   = 10                    % Frequency to calculate fields

ai = 22.86                      % Radius of input waveguide in mm
ao = 80%63.75                      % Radius of output waveguide in mm
bi = 10.16
bo = 60%49.78

pitch_fraction      = 10        % Choose a fraction between 10 to 5 (lambda_c / pitch_fraction)
delta               = 0.8       % Pitch to width ratio 0.7 to 0.9
wg_length           = 90;       % Length of feeding waveguide
num_of_corrugations = 40;
straight_width      = 2
cap_width           = 2

exc_mode = 'TE10';

SHOW_STRUCTURE_FIGURES  = 1;
RUN_SIMULATION          = 1;
PLOT_OUTPUT_SAME_WINDOW = 0;
USE_MODE_CONVERTER      = 0;
USE_WU_PROFILE          = 1;
USE_CORRUGATIONS        = 1;
SUBSTRACT_LEFTOVERS     = 1;

TIME_STEPS  = 10000
n_cell      = 50                    % cell size: lambda/n_cell
%%__________________________ END OF USER EDITABLE PARAMETERS __________________________

% Calculate center frequency fc based on narrow or wide bandwidth.
fratio = fmax/fmin             % ratio of fmax/fmin
if (fratio >= 2.4);            % check fmax/fmin is less than 2.4
    disp('Error, fmax/fmin is greater than 2.4!');
    fc = 0
    elseif (fratio <= 1.4);      % Use Narrowband formula if fmax <= 1.4fmin
        fc = sqrt(fmin*fmax)
        elseif (fmax >= 1.4*fmin && fmax <= 2.4*fmin); % Use wideband formula for 1.4fmin<=fmax<=2.4fmin
            fc = 1.2*fmin
endif

if (fratio <= 1.4);
    fo = 1.02*fc    % For narrow band choose fc <= fo <= 1.05fc
    elseif (fmax >= 1.4*fmin && fmax <= 2.4*fmin);
        fo = 1.10*fc    % For wideband choose 1.05fc <= fo <= 1.15fc
endif

unit = 1e-3;                    % Units in mm
lambda_c = 300/fcalc               % Center frequency wavelength
lambda_o = 300/fo               % Output frequency
p = lambda_c/pitch_fraction;    % Pitch in mm, lambda_c/10 to lambda_c/5
length = num_of_corrugations*p  % Length of horn profile
N = length/p                    % Total number of corrugations
z = 0:p:length;                 % z index distance array from 0 to length of horn


%%_____________________________ START OF 2D FIGURES DESIGN _____________________________

%%% Linear profile A %%%
a_profile = ai+(ao/2-ai/2)*z/length;

if (SHOW_STRUCTURE_FIGURES);
    figure
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
    d = 1:N+1;
    depth_step = (((300/fcalc)/2) - ((300/fcalc)/4)) / N

    for i = 1:N+1;
        d(i) = ((300/fcalc)/2) - i*depth_step;
    endfor
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
    z_for_a_profile(i+n+2) = z_for_a_profile(i+n)+delta*p;
    z_for_a_profile(i+n+3) = z_for_a_profile(i+n+2);
    z_for_a_profile(i+n+4) = z_for_a_profile(i+n+3)+(1-delta)*p;
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

% Add vertical surface at horn aperture
z_for_a_profile = [z_for_a_profile, z_for_a_profile(z_number)];
y_for_a_profile = [y_for_a_profile, a_offset(N)];
radmsh=y_for_a_profile;                 % radmesh to fix mesh lines to corrugations
% figure;                    % Uncomment these three lines for debugging
% plot(z_for_a_profile, y_for_a_profile);
% axis equal;

% Flip outer surface profile so that widest horn dimensions comes next in the outline coordinates
outer_surface = fliplr(a_offset);
z_flip = fliplr(z);
extent = z_for_a_profile(end);  % Fudge to make horn aperture planar for ring loaded slot MC
z_flip(1) = extent; % Fudge to make horn aperture planar for ring loaded slot MC
% Add outer profile and circular waveguide to horn
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
z_for_a_subs_profile = [z_flip,         -wg_length-2,       -wg_length-2,   0,  z_flip(1),  z_flip(1)];
y_for_a_subs_profile = [outer_surface,  ai+(lambda_c/2+2),  ao,             ao, ao,         outer_surface(1)];


%%% Linear profile B %%%
b_profile = bi+(bo/2-bi/2)*z/length;

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
    d = 1:N+1;
    depth_step = (((300/fcalc)/2) - ((300/fcalc)/4)) / N

    for i = 1:N+1;
        d(i) = ((300/fcalc)/2) - i*depth_step;
    endfor
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
    z_for_b_profile(i+n+2) = z_for_b_profile(i+n)+delta*p;
    z_for_b_profile(i+n+3) = z_for_b_profile(i+n+2);
    z_for_b_profile(i+n+4) = z_for_b_profile(i+n+3)+(1-delta)*p;
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
% a_offset is the inner horn profile shifted up to give the horn a thickness
a_offset = b_profile.+(lambda_c/2+2);
% figure;                    % Uncomment these three lines for debugging
% plot(z, a_offset);
% axis equal;

% Add vertical surface at horn aperture
z_for_b_profile = [z_for_b_profile, z_for_b_profile(z_number)];
y_for_b_profile = [y_for_b_profile, a_offset(N)];
radmsh=y_for_b_profile;                 % radmesh to fix mesh lines to corrugations
% figure;                    % Uncomment these three lines for debugging
% plot(z_for_b_profile, y_for_b_profile);
% axis equal;

% Flip outer surface profile so that widest horn dimensions comes next in the outline coordinates
outer_surface = fliplr(a_offset);
z_flip = fliplr(z);
extent = z_for_b_profile(end);  % Fudge to make horn aperture planar for ring loaded slot MC
z_flip(1) = extent; % Fudge to make horn aperture planar for ring loaded slot MC
% Add outer profile and circular waveguide to horn
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
endif

% Substraction volume for B
z_for_b_subs_profile = [z_flip,         -wg_length-2,       -wg_length-2,   0,  z_flip(1),  z_flip(1)];
y_for_b_subs_profile = [outer_surface,  bi+(lambda_c/2+2),  bo,             bo, bo,         outer_surface(1)];


%% Generate end cap to prevent the radiation coming out of the back of the horn
% Extract the end point of the structure
cap_point_z_1 = [z_flip, -wg_length](end);
cap_point_y_1 = [outer_surface, ai+(lambda_c/2+2)](end) - ai/2;
% From here we have to move in "y" to to met the "y" end of the cap.
cap_point_z_2 = cap_point_z_1;
cap_point_y_2 = -cap_point_y_1;
% We defined here the width of the cap
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
mesh.x = [(-a_offset(end)-(9*max_res)-lambda_max) -radmsh(1:4:end) 0 radmsh(1:4:end) (a_offset(end)+(9*max_res)+
                                                                                                        lambda_max)];
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.5); % Create a smooth mesh between specified fixed mesh lines
mesh.y = mesh.x;                                 % Same as x mesh
% Create fixed lines for the simulation box,port and given number of lines inside the horn
mesh.z = [-wg_length-lambda_max-(9*max_res) -wg_length-1 -wg_length -wg_length+10 0 z_for_a_profile(1:2:z_number) length+2*lambda_max+(9*max_res)];
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );

CSX = DefineRectGrid( CSX, unit, mesh );

%% ----->> Create Horn Geometry <<-----
%% Horn + waveguide
CSX = AddMetal(CSX, 'Corrugated_Horn');
CSX = AddMaterial( CSX, 'Air' );
CSX = SetMaterialProperty( CSX, 'Air', 'Epsilon', 1, 'Mue', 1 );

corrugated_coords_a = [y_for_a_profile; z_for_a_profile];
substract_coords_a  = [y_for_a_subs_profile; z_for_a_subs_profile];
corrugated_coords_b = [y_for_b_profile; z_for_b_profile];
substract_coords_b  = [y_for_b_subs_profile; z_for_b_subs_profile];
guard = 2

% Corrugated walls A
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 5, 1, 0, corrugated_coords_a, bo+(300/fcalc)+guard, 'Transform', {
    'Rotate_Y', -pi/2, 'Translate',[ num2str(ai/2) ',' num2str(-bo/2-(300/fcalc)/2+guard/2) ',0']
    });
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 5, 1, 0, corrugated_coords_a, bo+(300/fcalc)+guard, 'Transform', {
    'Rotate_Y', pi/2, 'Rotate_X', pi, 'Translate',[ num2str(-ai/2) ',' num2str(bo/2+(300/fcalc)/2+guard/2) ',0']
    });

if (SUBSTRACT_LEFTOVERS);
    CSX = AddLinPoly( CSX, 'Air', 10, 1, 0, substract_coords_a, bo+(300/fcalc)+guard, 'Transform', {
        'Rotate_Y', -pi/2, 'Translate',[ num2str(ai/2) ',' num2str(-bo/2-(300/fcalc)/2+guard/2) ',0']
        });
    CSX = AddLinPoly( CSX, 'Air', 10, 1, 0, substract_coords_a, bo+(300/fcalc)+guard, 'Transform', {
        'Rotate_Y', pi/2, 'Rotate_X', pi, 'Translate',[ num2str(-ai/2) ',' num2str(bo/2+(300/fcalc)/2+guard/2) ',0']
        });
endif

% Corrugated walls B
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 5, 1, 0, corrugated_coords_b, ao+(300/fcalc)+guard, 'Transform', {
    'Rotate_Y', -pi/2, 'Rotate_Z', -pi/2, 'Translate',[ num2str(-ao/2-(300/fcalc)/2+guard/2) ',' num2str(-bi/2) ',0']
    });
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 5, 1, 0, corrugated_coords_b, ao+(300/fcalc)+guard, 'Transform', {
    'Rotate_Y', -pi/2, 'Rotate_Z', pi/2, 'Translate',[ num2str(ao/2+(300/fcalc)/2+guard/2) ',' num2str(bi/2) ',0']
    });

if (SUBSTRACT_LEFTOVERS);
    CSX = AddLinPoly( CSX, 'Air', 10, 1, 0, substract_coords_b, ao+(300/fcalc)+guard, 'Transform', {
        'Rotate_Y', -pi/2, 'Rotate_Z', -pi/2, 'Translate',[ num2str(-ao/2-(300/fcalc)/2+guard/2) ',' num2str(-bi/2) ',0']
        });
    CSX = AddLinPoly( CSX, 'Air', 10, 1, 0, substract_coords_b, ao+(300/fcalc)+guard, 'Transform', {
        'Rotate_Y', -pi/2, 'Rotate_Z', pi/2, 'Translate',[ num2str(ao/2+(300/fcalc)/2+guard/2) ',' num2str(bi/2) ',0']
        });
endif


%% End cap to prevent the radiation coming out of the back of the horn
CSX = AddMetal(CSX, 'Cap');
cap_coords = [cap_y; cap_z];

CSX = AddLinPoly( CSX, 'Cap', 10, 1, 0, cap_coords, bo + 2 * straight_width, 'Transform',{
    'Rotate_Y', -pi/2, 'Translate',[ '0,' num2str(-bo/2 - straight_width) ',' num2str(-cap_width)]
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
Sim_Path = 'tmp';
Sim_CSX = 'Corrugated_Horn.xml';
[status, message, messageid] = rmdir( Sim_Path, 's'); % Clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % Create empty simulation folder

% Write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

% Show structure
CSXGeomPlot([Sim_Path '/' Sim_CSX], ['--export-STL=tmp']);

%% ----->> End of simulation folder <<----- 


%% ----->> Run openEMS <<-----
if(RUN_SIMULATION == 1)                                                                              %% Start Simulation
    %openEMS_opts = '--debug-PEC --no-simulation';   % Uncomment to visualise mesh in Paraview
    %RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts);
    RunOpenEMS(Sim_Path, Sim_CSX, '--numThreads=3');

    % Postprocessing & do the plots
    freq = linspace(f_start,f_stop,201);
    port = calcPort(port, Sim_Path, freq);

    Zin = port.uf.tot ./ port.if.tot;
    s11 = port.uf.ref ./ port.uf.inc;

    % Plot reflection coefficient S11
    figure
    if (PLOT_OUTPUT_SAME_WINDOW == 1);
        subplot (3, 2, 1)
    endif
    plot(freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2);
    xlim([fmin fmax]);
    ylim([-40 0]);
    set(gca, "linewidth",2, "fontsize", 14)
    grid on
    title('Reflection Coefficient S_{11}', 'FontSize', 16);
    xlabel('Frequency (GHz)','FontSize', 14);
    ylabel('Reflection Coefficient |S_{11}| (dB)','FontSize', 14);
    drawnow

    % NFFF plots

    % Calculate the far field at phi=0, 45 and at phi=90 degrees
    thetaRange = (0:0.2:359) - 180;
    disp('calculating far field at phi=[0 45 90] deg...');
    nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, [0 45 90]*pi/180);

    Dlog=10*log10(nf2ff.Dmax);      % Calculate maximum Directivity in dB

    % Display some antenna parameters from above calculations
    disp(['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
    disp(['directivity: Dmax = ' num2str(Dlog) ' dBi']);

    % Directivity
    if (PLOT_OUTPUT_SAME_WINDOW == 1);
        subplot (3, 2, 2)
    elseif
        figure
    endif
    plotFFdB(nf2ff,'xaxis','theta','param',[1 2 3]);
    ylim([-30 25]);
    xlim([-180 180]);
    grid on
    set(gca,"linewidth",2, "fontsize", 14, "XTick", -180:30:180, "YTick", -30:5:40)
    title(sprintf('Farfield Directivity @ %.2f GHz',fcalc),'FontSize', 16);
    xlabel('Theta (degrees)','FontSize', 14);
    ylabel('Directivity (dBi)','FontSize', 14);
    drawnow

    % Plot Ludwig3 cross polar
    if (PLOT_OUTPUT_SAME_WINDOW == 1);
        subplot (3, 2, 3)
    elseif
        figure
    endif
    plotFFcocx(nf2ff,'xaxis','theta','param',[2]);
    ylim([-30 25]);
    xlim([-180 180]);
    grid on
    set(gca,"linewidth",2, "fontsize", 14, "XTick", -180:30:180, "YTick", -30:5:40)
    title(sprintf('Farfield Directivity with Ludwig3 XPOL @ %.2f GHz',fcalc),'FontSize', 16);
    xlabel('Theta (degrees)','FontSize', 14);
    ylabel('Directivity (dBi)','FontSize', 14);
    drawnow

    % Polar plot
    if (PLOT_OUTPUT_SAME_WINDOW == 1);
        subplot (3, 2, 4)
    elseif
        figure
    endif
    leg=[];   %legend
    polarFF(nf2ff,'xaxis','theta','param',[1 2 3],'logscale',[-30 35], 'xtics', 12);
    title(sprintf('Farfield Directivity @ %.2f GHz',fcalc),'FontSize', 16);
    xlabel('Theta (degrees)','FontSize', 14);
    ylabel('Directivity (dBi)','FontSize', 14);
    drawnow

    %% Calculate 3D pattern
    %phiRange = sort(unique([-180:5:-100 -100:2.5:-50 -50:1:50 50:2.5:100 100:5:180]));
    %thetaRange = sort(unique([0:1:50 50:2:100 100:5:180]));
    phiRange = sort(unique([-180:1:-100 -100:1:-50 -50:1:50 50:1:100 100:1:180]));
    thetaRange = sort(unique([0:1:50 50:1:100 100:1:180]));

    disp('calculating 3D far field...');
    nf2ff = CalcNF2FF(nf2ff, Sim_Path, f0, thetaRange*pi/180, phiRange*pi/180, 'Verbose',2,'Outfile','nf2ff_3D.h5');

    if (PLOT_OUTPUT_SAME_WINDOW == 1);
        subplot (3, 2, 5)
    elseif
        figure
    endif
    colormap jet;
    plotFF3D(nf2ff, 'logscale', -40);        % plot 3D far field in dB

    % Save far field in VTK to plot in ParaView
    E_far_normalized = nf2ff.E_norm{1}/max(nf2ff.E_norm{1}(:));
    DumpFF2VTK([Sim_Path '/Farfield.vtk'],E_far_normalized,thetaRange,phiRange,'scale', 0.008, 'logscale', -30, 
                                                                                                    'maxgain', Dlog);
end                                                                                                    %% End Simulation
