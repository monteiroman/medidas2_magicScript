close all
clear
clc

%%_____________________________ USER EDITABLE PARAMETERS _____________________________
fmin    = 8                             % Minimum frequency in GHz
fmax    = 12                            % Maximum frequency in GHz
fcalc   = 10                            % Frequency to calculate fields
pitch_fraction = 8                      % Choose a fraction between 10 to 5 (lambda_c / pitch_fraction)
delta   = 0.8                           % Pitch to width ratio 0.7 to 0.9
sigma   = 0.42                          % Percentage factor for first slot depth, 0.4 to 0.5 
NMC     = 5                             % Number of corrugations in mode converter
wgl     = 30;                           % Length of circular feeding waveguide
num_of_corrugations = 60
corrugated_width    = 10.16
straight_width      = 2
cap_width           = 2


%                        aperture_wall_depth
%                               |   |
%                     _   __________    _
%                        |       ___|     aperture_wall_step
%                        |      |___    _
%  aperture_wall_height  |       ___|   _ aperture_wall_pitch = aperture_wall_pitch_fraction * aperture_wall_step
%                        |      |___
%                        |       ___|
%                        |      |
%                     _  |      |
%
%                        |      |
%                  aperture_wall_width  
%
aperture_wall_pitch_fraction = 0.5                                          % 0 to 1
aperture_wall_N_corrugations = 10
aperture_wall_step = 5                                                     % in mm
aperture_wall_pitch = aperture_wall_step*aperture_wall_pitch_fraction       
aperture_wall_height = aperture_wall_N_corrugations*aperture_wall_step      % Length of horn aperture walls
aperture_wall_width = 5
aperture_wall_depth = 5

exc_mode = 'TE10';

SHOW_STRUCTURE_FIGURES  = 1;
RUN_SIMULATION          = 1;
USE_MODE_CONVERTER      = 0;

TIME_STEPS = 50000
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
lambda_c = 300/fc               % Center frequency wavelength
lambda_o = 300/fo               % Output frequency
ai = 22.86%(3 * lambda_c)/(2*pi)      % Radius of input waveguide in mm
ao = 3.5*lambda_c%1.95*lambda_c              % Radius of output waveguide in mm
p = lambda_c/pitch_fraction;    % Pitch in mm, lambda_c/10 to lambda_c/5
length = num_of_corrugations*p  % Length of horn profile
N = length/p                    % Total number of corrugations
kc = (2*pi)/lambda_c            % Wave number at center frequency
ko = (2*pi)/lambda_o            % Wave number at output frequency
z = 0:p:length;                 % z index distance array from 0 to length of horn
r_app = ao*1e-3;                % Aperture radius
A_app = pi*(r_app)^2;           % Aperture area for gain calculation

aperture_wall_length = ao*2;

%%_____________________________ START OF 2D FIGURES DESIGN _____________________________
%%% Linear profile %%%
a = ai+(ao-ai)*z/length;

if (SHOW_STRUCTURE_FIGURES);
    figure
    subplot (3, 2, 1)
    plot(z, a);    
    set(gca, "linewidth",2, "fontsize", 14 )
    xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
    ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
    title( 'Linear Horn Profile', 'FontSize', 16 );
    axis equal;
endif

if (USE_MODE_CONVERTER == 1);
    % Mode Converter depths for element j
    ajmc = a(1:NMC);                     % Index range for mode converter
    idx = 1:NMC;
    djmc = (sigma-((idx-1)./NMC).*(sigma-(0.25.*exp(1./(2.114.*(kc*ajmc).^1.134)))))*lambda_c;
    % Depth of remaining corrugations
    aj = a(NMC+1:end);
    idx = NMC+1:N+1;
    dj = ((lambda_c/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-NMC-1)/(N-NMC-1)).*((lambda_c/4).*exp(
                                            1./(2.114.*(kc*aj).^1.134))-(lambda_o/4).*exp(1./(2.114.*(ko*ao).^1.134)));
    d = [djmc, dj];       % Combining the mode converter and horn depth values
else
    % Depth of remaining corrugations
    aj = a;
    idx = 1:N+1;
    dj = ((lambda_c/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-1)/(N-1)).*((lambda_c/4).*exp(
                                            1./(2.114.*(kc*aj).^1.134))-(lambda_o/4).*exp(1./(2.114.*(ko*ao).^1.134)));
    d = dj;       % Combining the mode converter and horn depth values
endif

% Generate z,y coordinates as len and rad vector
n = 0;
len(1) = 0;
len(2) = 0;
for i = 1:N;
    rad(i+n) = a(i);
    rad(i+n+1) = a(i)+d(i);
    rad(i+n+2) = a(i)+d(i);
    rad(i+n+3) = a(i+1);
    rad(i+n+4) = a(i+1);
    len(i+n+2) = len(i+n)+delta*p;
    len(i+n+3) = len(i+n+2);
    len(i+n+4) = len(i+n+3)+(1-delta)*p;
    len(i+n+5) = len(i+n+4);
    n = n+3;
endfor

z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
len = len(1:z_number); % Truncate z axis data points to equal rad vector length

if (SHOW_STRUCTURE_FIGURES);
    subplot (3, 2, 2)
    plot(len,rad);
    set(gca, "linewidth",2, "fontsize", 14 )
    xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
    ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
    title( 'Mode Converter and Corrugation Profile', 'FontSize', 16 );
    axis equal;   % Scale axis equally for aspect ratio 1:1
endif

% Add the rest of the geometry to create a closed path
% a_offset is the inner horn profile shifted up to give the horn a thickness
a_offset = a.+(lambda_c/2+2);
% figure;                    % Uncomment these three lines for debugging
% plot(z, a_offset);
% axis equal;

% Add vertical surface at horn aperture
len = [len, len(z_number)];
rad = [rad, a_offset(N)];
radmsh=rad;                 % radmesh to fix mesh lines to corrugations
% figure;                    % Uncomment these three lines for debugging
% plot(len, rad);
% axis equal;

% Flip outer surface profile so that widest horn dimensions comes next in the outline coordinates
outer_surface = fliplr(a_offset);
z_flip = fliplr(z);
extent = len(end);  % Fudge to make horn aperture planar for ring loaded slot MC
z_flip(1) = extent; % Fudge to make horn aperture planar for ring loaded slot MC
% Add outer profile and circular waveguide to horn
len = [len, z_flip,         -wgl,               -wgl,   0];
rad = [rad, outer_surface,  ai+(lambda_c/2+2),  ai,     ai];

if (SHOW_STRUCTURE_FIGURES);
    subplot (3, 2, 3)
    plot(len,rad);
    set(gca, "linewidth",2, "fontsize", 14 )
    xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
    ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
    title( 'Complete Corrugated Horn Profile', 'FontSize', 16 );
    axis equal;   % Scale axis equally for aspect ratio 1:1
endif

% Determine straight faces
straight_len = [z_flip, -wgl, -wgl, 0, length, length, length];
straight_rad = [outer_surface, ai+(lambda_c/2+2), ai/2, ai/2, ai/2, ao-ai/2, outer_surface(1)];

if (SHOW_STRUCTURE_FIGURES);
    subplot (3, 2, 4)
    plot(straight_len, straight_rad);
    set(gca, "linewidth",2, "fontsize", 14 )
    xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
    ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
    title( 'Straight walls', 'FontSize', 16 );
endif

% Determine end cap to prevent the radiation coming out of the back of the horn
% Extract the end point of the structure
cap_point_z_1 = [z_flip, -wgl](end);
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

if (SHOW_STRUCTURE_FIGURES);
    subplot (3, 2, 5)
    plot(cap_z, cap_y, 'o-r');
    set(gca, "linewidth",2, "fontsize", 14 )
    xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
    ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
    title( 'Cap', 'FontSize', 16 );
    axis equal;   % Scale axis equally for aspect ratio 1:1
endif

% Define the apertue walls
aperture_wall_x(1) = 0;
aperture_wall_y(1) = 0;
aperture_wall_x(2) = 0;
aperture_wall_y(2) = aperture_wall_height;
aperture_wall_x(3) = aperture_wall_width;
aperture_wall_y(3) = aperture_wall_y(2);

n = 0;
for i = 3:aperture_wall_N_corrugations+2;
    aperture_wall_x(i+n+1) = aperture_wall_x(i+n) + aperture_wall_depth;
    aperture_wall_y(i+n+1) = aperture_wall_y(i+n);
    aperture_wall_x(i+n+2) = aperture_wall_x(i+n+1);
    aperture_wall_y(i+n+2) = aperture_wall_y(i+n+1) - aperture_wall_pitch;
    aperture_wall_x(i+n+3) = aperture_wall_x(i+n+2) - aperture_wall_depth;
    aperture_wall_y(i+n+3) = aperture_wall_y(i+n+2);
    aperture_wall_x(i+n+4) = aperture_wall_x(i+n+3);
    aperture_wall_y(i+n+4) = aperture_wall_y(i+n+3) - (aperture_wall_step - aperture_wall_pitch);
    n = n + 3;
endfor

% Set corrugations to zero in x and y.
aperture_wall_x = aperture_wall_x - (aperture_wall_width + aperture_wall_depth);
aperture_wall_y = aperture_wall_y - (aperture_wall_width + aperture_wall_depth);


if (SHOW_STRUCTURE_FIGURES);
    subplot (3, 2, 6)
    plot(aperture_wall_x, aperture_wall_y);    
    set(gca, "linewidth",2, "fontsize", 14 )
    xlabel( 'Dimension in x Direction (mm)', 'FontSize', 14 );
    ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
    title( 'Aperture walls', 'FontSize', 16 );
    axis equal;   % Scale axis equally for aspect ratio 1:1
endif


%%%_____________________________ END OF 2D FIGURES DESIGN _____________________________
%

% openEMS setup begins here
% EM related physical constants
physical_constants;

% frequency range of interest
f_start =  fmin*1e9;
f_stop  =  fmax*1e9;

% frequency to calculate fields
f0 = 12.46*1e9;

%% setup FDTD parameter & excitation function
FDTD = InitFDTD( 'NrTS', 50000, 'EndCriteria', 0.5e-3 );
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond(FDTD, BC);

%% setup CSXCAD geometry & mesh
max_res = c0/(f_stop)/unit/20; % cell size: lambda/20
CSX = InitCSX();               % Initialise CSX structure

% Calculate lambda/4 at lowest frequency to use as distance to nf2ff surfaces
lambda_max = c0/f_start/unit/4;

% Create fixed lines for the simulation box, structure and port
mesh.x = [(-a_offset(end)-(9*max_res)-lambda_max) -radmsh(1:4:end) 0 radmsh(1:4:end) (a_offset(end)+(9*max_res)+
                                                                                                        lambda_max)];
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.5); % Create a smooth mesh between specified fixed mesh lines
mesh.y = mesh.x/2;                                 % Same as x mesh
% Create fixed lines for the simulation box,port and given number of lines inside the horn
mesh.z = [-wgl-lambda_max-(9*max_res) -wgl-1 -wgl -wgl+10 0 len(1:2:z_number) length+2*lambda_max+(9*max_res)];
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );

CSX = DefineRectGrid( CSX, unit, mesh );

%% ----->> Create Horn Geometry <<-----

%% Horn + waveguide
CSX = AddMetal(CSX, 'Corrugated_Horn');
corrugated_coords = [rad; len];
straight_coords = [straight_rad; straight_len];
aperture_wall_coords = [aperture_wall_x; aperture_wall_y];

% Corrugated walls
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 10, 1, 0, corrugated_coords, corrugated_width, 'Transform', {
    'Rotate_Y', -pi/2, 'Translate',[ num2str(ai/2) ',' num2str(-corrugated_width/2) ',0']
    });
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 10, 1, 0, corrugated_coords, corrugated_width, 'Transform', {
    'Rotate_Y', pi/2, 'Rotate_X', pi, 'Translate',[ num2str(-ai/2) ',' num2str(corrugated_width/2) ',0']
    });
% Straight walls
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 10, 1, 0, straight_coords, straight_width, 'Transform', {
    'Rotate_Y', -pi/2, 'Translate',[ num2str(ai/2) ',' num2str(-corrugated_width/2-straight_width) ',0']
    });
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 10, 1, 0, straight_coords, -straight_width, 'Transform', {
    'Rotate_Y', pi/2, 'Rotate_X', pi, 'Translate',[ num2str(-ai/2) ',' num2str(-corrugated_width/2-straight_width) ',0']
    });

CSX = AddLinPoly( CSX, 'Corrugated_Horn', 10, 1, 0, straight_coords, -straight_width, 'Transform', {
    'Rotate_Y', -pi/2, 'Translate',[ num2str(ai/2) ',' num2str(corrugated_width/2+straight_width) ',0']
    });
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 10, 1, 0, straight_coords, straight_width, 'Transform', {
    'Rotate_Y', pi/2, 'Rotate_X', pi,  'Translate',[ num2str(-ai/2) ',' num2str(corrugated_width/2+straight_width) ',0']
    });
% Corrugated aperture walls
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 10, 1, 0, aperture_wall_coords, aperture_wall_length, 'Transform',{
    'Rotate_Z', pi/2,'Translate',[num2str(ao) ',' num2str(corrugated_width+straight_width*2) ',' num2str(length)]
    });
CSX = AddLinPoly( CSX, 'Corrugated_Horn', 10, 1, 0, aperture_wall_coords, aperture_wall_length, 'Transform',{
    'Rotate_Z', -pi/2,'Translate',[num2str(-ao) ',' num2str(-corrugated_width-straight_width*2) ',' num2str(length)]
    });

%% End cap to prevent the radiation coming out of the back of the horn
CSX = AddMetal(CSX, 'Cap');
cap_coords = [cap_y; cap_z];

CSX = AddLinPoly( CSX, 'Cap', 10, 1, 0, cap_coords, corrugated_width + 2 * straight_width, 'Transform',{
    'Rotate_Y', -pi/2, 'Translate',[ '0,' num2str(-corrugated_width/2 - straight_width) ',' num2str(-cap_width)]
    });
%% ----->> End of model geometry <<-----

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
                                                                                                        % End simulate
