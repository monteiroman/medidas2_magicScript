%
% EXAMPLE / waveguide / Coax to Waveguide Adapter
%
% This example demonstrates:
%  - mixed port types

% Tested with
%  - Octave 4.0.0
%  - openEMS v0.0.32
%
close all
clear
clc

%% switches
postproc_only = 0;

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 2.54/100; %drawing unit in inches
NumTS = 50000; %max. number of timesteps

% waveguide dimensions and mode
m = 1;
n = 0;

length = 1.0;
a = .75;     %waveguide width  WR-75 works
b = .375;    %waveguide heigth
WallThickness = 0.04;   % thickness of top wall
BackShort = .235;    % distance from short to center of probe

% SMA dimensions
InnerCondSMA = 0.05;      %inner diameter (inch)
OuterCondSMA = 0.166;     %inner diam of outer conductor (inch)
OuterCondODSMA = 0.190;    % outer diam of outer conductor (inch)
ProbeDepth = 0.195;       % Probe insertion depth inside waveguide
SMALength = 0.4;          % length of SMA connector (inch)
epsR = 2.08;

mesh_res = [.007 .01 .007];

f_start = 10e9;    
f_stop  = 15e9;   
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

if (postproc_only==0)
    [status, message, messageid] = rmdir(Sim_Path,'s');
    [status, message, messageid] = mkdir(Sim_Path);
end

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD('NrTS',NumTS,'EndCriteria', 1e-5); 
FDTD = SetGaussExcite(FDTD,.5*(f_stop+f_start),.5*(f_stop-f_start));    

BC = {'PEC' 'PEC' 'PEC' 'PML_8' 'PML_8' 'PEC'};   
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CSX = InitCSX();

mesh.x = SmoothMeshLines2( [a/2-OuterCondODSMA/2 a/2-InnerCondSMA/2 a/2+InnerCondSMA/2 a/2+OuterCondODSMA/2],mesh_res(1)/3,1.3);
mesh.x = SmoothMeshLines([0 mesh.x a], mesh_res(1),1.3);

mesh.y = SmoothMeshLines2([b-ProbeDepth b b+WallThickness],mesh_res(2)/2,1.3);
mesh.y = SmoothMeshLines([0 mesh.y b+WallThickness+SMALength], mesh_res(2));

mesh.z = SmoothMeshLines2([length-BackShort-OuterCondODSMA/2 length-BackShort-InnerCondSMA/2 length-BackShort+InnerCondSMA/2,length-BackShort+OuterCondODSMA/2],mesh_res(3)/3,1.3);
mesh.z = SmoothMeshLines([0 mesh.z length], mesh_res(3));

CSX = DefineRectGrid(CSX, unit,mesh);

%% add metal top to waveguide   
CSX = AddMetal(CSX,'metal');  %% metal is PEC
CSX = AddMetal(CSX,'metal2');  %% metal is PEC

CSX = AddBox(CSX,'metal',0, [0 b 0], [a b+WallThickness length]);

%% drill hole for SMA probe to enter
CSX = AddMaterial(CSX,'air');
CSX = SetMaterialProperty(CSX,'air','Epsilon',1.0,'Mue',1.0);
CSX = AddCylinder(CSX,'air',20,[a/2 b length-BackShort], [a/2 b+WallThickness length-BackShort],OuterCondSMA/2);

%% add SMA center conductor probe to waveguide
CSX = AddCylinder(CSX,'metal2',30, [a/2 b-ProbeDepth length-BackShort], [a/2 b+WallThickness length-BackShort],InnerCondSMA/2);

%% add teflon dielectric material

CSX = AddMaterial( CSX, 'teflon' );
CSX = SetMaterialProperty( CSX, 'teflon', 'Epsilon', epsR);
%% no need for AddBox because teflon shape is set in AddCoaxialPort

%%% coax and coax port #1

start = [a/2,b+WallThickness+SMALength,length-BackShort];
stop  = [a/2,b+WallThickness,length-BackShort];     
[CSX,port{1}] = AddCoaxialPort( CSX, 0, 1, 'metal', 'teflon', start, stop, 'y',InnerCondSMA/2 , OuterCondSMA/2, OuterCondODSMA/2,'ExciteAmp', 1,'FeedShift', 10*mesh_res(2) );

% waveguide port 2

start = [mesh.x(1)   mesh.y(1)   mesh.z(14)];
stop  = [mesh.x(end) b mesh.z(15)];
[CSX, port{2}] = AddRectWaveGuidePort( CSX, 0, 2, start, stop, 'z', a*unit, b*unit, TE_mode);

%% define dump box... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,4');
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

%% Write openEMS compatible xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (postproc_only==0)
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);
    CSXGeomPlot([Sim_Path '/' Sim_CSX], ['--export-STL=tmp']);
 
    RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts, Settings)
end

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
figure
dump_file = [Sim_Path '/Et.h5'];
PlotArgs.slice = {a/2*unit (b/2)*unit (.05)*unit};
PlotArgs.pauseTime=0.01;
PlotArgs.component=0;
PlotArgs.Limit = 'auto';
PlotHDF5FieldData(dump_file, PlotArgs);

disp('done...')
