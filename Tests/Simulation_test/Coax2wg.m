%
% Adaptador Conector N a Guía de onda WR-90 -> Keysight/Agilent x281A
%
close all
clear
clc

%% switches
postproc_only = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulación %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
physical_constants;
unit = 1e-3; %unidades en milímetros
NumTS = 30000; %max. number of timesteps

% Modo de transmisión TE10
m = 1;
n = 0;

clearance=6;
cl_x=10;
cl_y=4;

%%%%% Dimensiones de la guía %%%%%
length = 19.8;
b = 22.86;     %waveguide width  WR-75 works
a = 10.16;    %waveguide heigth
WallThickness = 5.48;   % thickness of top wall
BackShort = 4.3;    % distance from short to center of probe


%%%%% Dimensiones de la tapa %%%%%%
LadoTapa =41.4;
DepthTapa=5.6;


%%%%% Conector N %%%%%
InnerCondN = 3.086;      %inner diameter
OuterCondN = 8.636;     %inner diam of outer conductor
OuterCondODN = 16;    % outer diam of outer conductor
ProbeDepth = a/2;       % Probe insertion depth inside waveguide  ?????
NLength = 10.72;          % length of N connector
epsR = 2.08;


%%%%% Setup %%%%%%

f_start = 9e9;    
f_stop  = 12e9;   
TE_mode = 'TE10';

mesh_res = [.3 .3 .3];
%mesh_res = (c0/f_stop)./[3 3 3];


freq = linspace(f_start,f_stop,101);

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

BC = {'PML_10' 'PML_10' 'PML_10' 'PML_10' 'PML_10' 'PML_10'}; % BC changed
FDTD = SetBoundaryCond(FDTD,BC);

%% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CSX = InitCSX();

v_x=clearance+cl_x+[0 WallThickness WallThickness+a/2-OuterCondODN/2 WallThickness+a/2-InnerCondN/2 WallThickness+a/2+InnerCondN/2 WallThickness+a/2+OuterCondODN/2 a+WallThickness a+2*WallThickness ];
mesh.x = SmoothMeshLines2(v_x,mesh_res(1)/3,1.3);
mesh.x = SmoothMeshLines([0 mesh.x clearance+cl_x+WallThickness+a/2 a+2*WallThickness+2*clearance+cl_x LadoTapa+2*WallThickness], mesh_res(1),1.3);

v_y=clearance+cl_y+[0 WallThickness b+WallThickness b+2*WallThickness];
mesh.y = SmoothMeshLines2(v_y,mesh_res(2)/2,1.3);
mesh.y = SmoothMeshLines([0 mesh.y  b+2*WallThickness+2*clearance+cl_y b+2*WallThickness+3*clearance+cl_y], mesh_res(2));

v_z=clearance+[0 WallThickness length/2-OuterCondODN/2+WallThickness length/2-InnerCondN/2+WallThickness length/2+WallThickness length/2+InnerCondN/2+WallThickness length/2+OuterCondODN/2+WallThickness length+WallThickness];
mesh.z = SmoothMeshLines2(v_z,mesh_res(3)/3,1.3);
mesh.z = SmoothMeshLines([0 mesh.z clearance+length/2+WallThickness length+WallThickness+2*clearance length+2*WallThickness+2*clearance], mesh_res(3),1.3);

CSX = DefineRectGrid(CSX, unit,mesh);


%%%%%   Construyo la guía de onda   %%%%%
CSX = AddMetal(CSX,'metal');  %% metal is PEC
%top
CSX = AddBox(CSX,'metal',0, [clearance+WallThickness b+WallThickness+clearance clearance+WallThickness], [clearance+a+2*WallThickness b+2*WallThickness+clearance length+clearance+WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});
%bottom
CSX = AddBox(CSX,'metal',0, [clearance+WallThickness WallThickness+clearance clearance+WallThickness], [clearance+a+2*WallThickness clearance length+clearance+WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});
%left
CSX = AddBox(CSX,'metal',0, [clearance clearance clearance+WallThickness], [clearance+WallThickness b+2*WallThickness+clearance length+clearance+WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});
%right
CSX = AddBox(CSX,'metal',0, [clearance+a+WallThickness clearance clearance+WallThickness], [clearance+a+2*WallThickness b+2*WallThickness+clearance length+clearance+WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});
%back
CSX = AddBox(CSX,'metal',0, clearance*[1 1 1], clearance+[a+2*WallThickness b+2*WallThickness WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});

%%%%%   Agrego el conector N    %%%%%
% drill hole for N probe to enter
CSX = AddMaterial(CSX,'air');
CSX = SetMaterialProperty(CSX,'air','Epsilon',1.0,'Mue',1.0);
CSX = AddCylinder(CSX,'air',20,[a/2+WallThickness b+WallThickness WallThickness+BackShort], [a/2+WallThickness b+2*WallThickness+NLength WallThickness+BackShort],OuterCondN/2,'Transform',{'Translate', [clearance+cl_x,clearance+cl_y,clearance]});
% add N center conductor probe to waveguide
CSX = AddCylinder(CSX,'metal',30, [a/2+WallThickness b+WallThickness-ProbeDepth WallThickness+BackShort], [a/2+WallThickness b+2*WallThickness WallThickness+BackShort],InnerCondN/2,'Transform',{'Translate', clearance+[cl_x,cl_y,0]});

% add teflon dielectric material

CSX = AddMaterial( CSX, 'teflon' );
CSX = SetMaterialProperty( CSX, 'teflon', 'Epsilon', epsR);
% no need for AddBox because teflon shape is set in AddCoaxialPort

%%%%%   Construyo la tapa   %%%%%
%top
CSX = AddBox(CSX,'metal',0, [clearance+WallThickness+a/2 clearance+WallThickness+b clearance+length+WallThickness], [clearance+WallThickness+a/2+LadoTapa/2 clearance+WallThickness+b/2+LadoTapa/2 clearance+length+2*WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});
CSX = AddBox(CSX,'metal',0, [clearance+WallThickness+a/2 clearance+WallThickness+b clearance+length+WallThickness], [clearance+WallThickness+a/2-LadoTapa/2 clearance+WallThickness+b/2+LadoTapa/2 clearance+length+2*WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});
% %bottom
CSX = AddBox(CSX,'metal',0, [clearance+WallThickness+a/2 clearance+WallThickness clearance+length+WallThickness], [clearance+WallThickness+a/2+LadoTapa/2 clearance+WallThickness+b/2-LadoTapa/2 clearance+length+2*WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});
CSX = AddBox(CSX,'metal',0, [clearance+WallThickness+a/2 clearance+WallThickness clearance+length+WallThickness], [clearance+WallThickness+a/2-LadoTapa/2 clearance+WallThickness+b/2-LadoTapa/2 clearance+length+2*WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});
%left
CSX = AddBox(CSX,'metal',0, [clearance+WallThickness clearance+WallThickness clearance+WallThickness+length], [clearance+WallThickness+a/2-LadoTapa/2 b+WallThickness+clearance length+clearance+2*WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});
%right
CSX = AddBox(CSX,'metal',0, [clearance+WallThickness+a clearance+WallThickness clearance+length+WallThickness], [clearance+WallThickness+a/2+LadoTapa/2 clearance+WallThickness+b clearance+length+2*WallThickness],'Transform',{'Translate', [cl_x,cl_y,0]});

%coax port #1

start = clearance+[cl_x+a/2+WallThickness,cl_y+b+2*WallThickness+NLength,length-BackShort-WallThickness];
stop  = clearance+[cl_x+a/2+WallThickness,cl_y+b+1*WallThickness,length-BackShort-WallThickness];     
[CSX,port{1}] = AddCoaxialPort( CSX, 0, 1, 'metal', 'teflon', start, stop, 'y',InnerCondN/2 , OuterCondN/2, OuterCondODN/2,'ExciteAmp', 1,'FeedShift', 10*mesh_res(2));

% waveguide port 2

start = [clearance+cl_x+WallThickness   clearance+cl_y+WallThickness   mesh.z(297)];
stop  = [clearance+cl_x+WallThickness+a clearance+cl_y+WallThickness+b mesh.z(300)];
[CSX, port{2}] = AddRectWaveGuidePort( CSX, 0, 2, start, stop, 'z', a*unit, b*unit, TE_mode);

%% define dump box... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,4');
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

%% Write openEMS compatible xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (postproc_only==0)
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);
    CSXGeomPlot([Sim_Path '/' Sim_CSX]);
 
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