close all
clear
clc

addpath(genpath("../../src"));

YES = 1;
NO  = 0;
ON  = 1;
OFF = 0;

%% switches
run_simulation = YES;

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit            = 1e-3; 
NumTS = 50000; %max. number of timesteps

% waveguide dimensions and mode
m = 1;
n = 0;

length          = 60;
b               = 10.16;            % waveguide width  WR-75 works
a               = 22.86;            % waveguide heigth
WallThickness   = 2;               % thickness of top wall
BackShort       = 6;                % distance from short to center of probe


%%%%% Conector N %%%%%
InnerCond_N     = 3.086;            %inner diameter
OuterCond_N     = 8.636;            %inner diam of outer conductor
OuterCondOD_N   = 14;               % outer diam of outer conductor
ProbeDepth      = 5.5;              % Probe insertion depth inside waveguide  ?????
N_Length        = 10.72;            % length of N connector
epsR            = 2.08;

space = 5;

mesh_res = [.5 .5 .5];

%%%%% Setup %%%%%%

f_start = 8e9;    
f_stop  = 12e9;   
TE_mode = 'TE10';

freq = linspace(f_start,f_stop,201);

k = 2*pi*freq/c0;
kc = sqrt((m*pi/a/unit)^2 + (n*pi/b/unit)^2);
fc = c0*kc/2/pi;          %cut-off frequency
beta = sqrt(k.^2 - kc^2); %waveguide phase-constant
ZL_a = k * Z0 ./ beta;    %analytic waveguide impedance

disp([' Cutoff frequencies for this mode and waveguide is: ' num2str(fc/1e9) ' GHz']);

if (f_start<fc)
    warning('openEMS:example','f_start is smaller than the cutoff-frequency, this may result in a long simulation... ');
end


%% define and openEMS options for RunOpenEMS.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '';
% openEMS_opts = [openEMS_opts ' --disable-dumps'];
% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --engine=basic'];

Settings = [];
Settings.LogFile = 'openEMS.log';

Sim_Path = 'tmp';
Sim_CSX = 'coax2wg.xml';

confirm_recursive_rmdir(0);                                 % Do not asks if remove directory
[status, message, messageid] = rmdir(Sim_Path,'s');
[status, message, messageid] = mkdir(Sim_Path);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD('NrTS',NumTS,'EndCriteria', 1e-5); 
FDTD = SetGaussExcite(FDTD,.5*(f_stop+f_start),.5*(f_stop-f_start));    

BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};   
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CSX = InitCSX();

mesh.x = SmoothMeshLines2([0 a/2-OuterCondOD_N/2 a/2-InnerCond_N/2 a/2+InnerCond_N/2 a/2+OuterCondOD_N/2 a],mesh_res(1)/3,1.3);
mesh.x = SmoothMeshLines([-(space+WallThickness) mesh.x a+(space+WallThickness)], mesh_res(1),1.3);

mesh.y = SmoothMeshLines2([b-ProbeDepth b -b+WallThickness],mesh_res(2)/2,1.3);
mesh.y = SmoothMeshLines([-(space+WallThickness) mesh.y b+WallThickness+N_Length+(space+WallThickness)], mesh_res(2));

mesh.z = SmoothMeshLines2([length-BackShort-OuterCondOD_N/2 length-BackShort-InnerCond_N/2 length-BackShort+InnerCond_N/2,length-BackShort+OuterCondOD_N/2],mesh_res(3)/3,1.3);
mesh.z = SmoothMeshLines([0 mesh.z length+(space+WallThickness)], mesh_res(3));

CSX = DefineRectGrid(CSX, unit,mesh);

%% add metal top to waveguide   
CSX = AddMetal(CSX,'metal');  %% metal is PEC
CSX = AddMetal(CSX,'metal2');  %% metal is PEC
% AddBox info at:
%       https://openems.de/index.php/Box.html
CSX = AddBox(CSX,'metal',0, [0 b 0], [a b+WallThickness length]);
CSX = AddBox(CSX,'metal',0, [0 0 0], [a -WallThickness length]);

CSX = AddBox(CSX,'metal',0, [a -WallThickness 0], [a+WallThickness b+WallThickness length]);
CSX = AddBox(CSX,'metal',0, [0 -WallThickness 0], [-WallThickness b+WallThickness length]);

CSX = AddBox(CSX,'metal',0, [-WallThickness -WallThickness length], [a+WallThickness b+WallThickness length+WallThickness]);


%% drill hole for _N probe to enter
CSX = AddMaterial(CSX,'air');
CSX = SetMaterialProperty(CSX,'air','Epsilon',1.0,'Mue',1.0);
% More info of AddCylinder here:
%                       https://openems.de/index.php/Cylinder.html
CSX = AddCylinder(CSX,'air',20,[a/2 b length-BackShort], [a/2 b+WallThickness length-BackShort],OuterCond_N/2);

%% add _N center conductor probe to waveguide
CSX = AddCylinder(CSX,'metal2',30, [a/2 b-ProbeDepth length-BackShort], [a/2 b+WallThickness length-BackShort],InnerCond_N/2);

%% add teflon dielectric material

CSX = AddMaterial( CSX, 'teflon' );
CSX = SetMaterialProperty( CSX, 'teflon', 'Epsilon', epsR);
%% no need for AddBox because teflon shape is set in AddCoaxialPort

%%% coax and coax port #1

start = [a/2,b+WallThickness+N_Length,length-BackShort];
stop  = [a/2,b+WallThickness,length-BackShort];     
[CSX,port{1}] = AddCoaxialPort( CSX, 0, 1, 'metal', 'teflon', start, stop, 'y',InnerCond_N/2 , OuterCond_N/2, OuterCondOD_N/2,'ExciteAmp', 1,'FeedShift', 10*mesh_res(2) );

% waveguide port 2

% start = [mesh.x(1)   mesh.y(1)   mesh.z(14)];
start = [0 0 mesh.z(14)];
% stop  = [mesh.x(end) b mesh.z(15)];
stop  = [a b mesh.z(15)];
[CSX, port{2}] = AddRectWaveGuidePort( CSX, 0, 2, start, stop, 'z', a*unit, b*unit, TE_mode);

%% define dump box... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,4');
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

%% Write openEMS compatible xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);
CSXGeomPlot([Sim_Path '/' Sim_CSX], ['--export-STL=tmp']);

if (run_simulation == YES)
    RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts, Settings)

    port = calcPort( port, Sim_Path, freq);

    % must correct s21 by ratio of each port impedance
    % impedance of the coax is still not perfect due to mesh discretization 
    % in a cartesian grid

    s11 = port{1}.uf.ref./ port{1}.uf.inc;
    s21 = sqrt(real(port{1}.ZL_ref)./port{2}.ZL_ref).*port{2}.uf.ref./port{1}.uf.inc;

    %% plot s-parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    plot(freq*1e-9,20*log10(abs(s11)),'k-','Linewidth',2);
    xlim([freq(1) freq(end)]*1e-9);
    grid on;
    hold on;
    plot(freq*1e-9,20*log10(abs(s21)),'r--','Linewidth',2);
    l = legend('S_{11}','S_{21}','Location','Best');
    set(l,'FontSize',12);
    ylabel('S-Parameter (dB)','FontSize',12);
    xlabel('frequency (GHz) \rightarrow','FontSize',12);

    %% plot impedance of coax port
    figure
    plot(freq*1e-9,real(port{1}.ZL_ref),'k-','Linewidth',2);
    xlim([freq(1) freq(end)]*1e-9);
    grid on;
    hold on;
    set(l,'FontSize',12);
    ylabel('coax impedance','FontSize',12);
    xlabel('frequency (GHz) \rightarrow','FontSize',12);


    %% Plot the field dumps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Great diagnostic tool
    % figure
    % dump_file = [Sim_Path '/Et.h5'];
    % PlotArgs.slice = {a/2*unit (b/2)*unit (.05)*unit};
    % PlotArgs.pauseTime=0.01;
    % PlotArgs.component=0;
    % PlotArgs.Limit = 'auto';
    % PlotHDF5FieldData(dump_file, PlotArgs);
end

disp('done...')
