close all
clear
clc

make_horn();

%% ----->> Run openEMS <<-----
 openEMS_opts = '--debug-PEC --no-simulation';   % Uncomment to visualise mesh in Paraview
 RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts);
%RunOpenEMS(Sim_Path, Sim_CSX);%, '--numThreads=3');

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
