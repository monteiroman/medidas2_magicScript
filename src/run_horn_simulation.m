%%
%%  Simulation run function.
%%
%%  Author: Tiago Monteiro
%%
function run_horn_simulation(Sim, port, nf2ff)
    %%
    %%  Runs the OpenEMS simulation and stores output graphics.
    %%
    %%  Parameters:
    %%          Sim: Horn simulation parameters.
    %%          port: Waveguide port for later simulation.
    %%          nf2ff: Near field to far field simulation box.
    %%
    %%  Returns:
    %%          void
    %%

    %% ----->> Run openEMS <<-----
    % openEMS_opts = '--debug-PEC --no-simulation';   % Uncomment to visualise mesh in Paraview
    % RunOpenEMS(Sim.Sim_Path, Sim.Sim_CSX, openEMS_opts);
    RunOpenEMS(Sim.Sim_Path, Sim.Sim_CSX, '--dump-statistics');%, '--numThreads=3');

    % frequency range of interest
    f_start =  Sim.fmin*1e9;
    f_stop  =  Sim.fmax*1e9;

    % frequency to calculate fields
    f0 = Sim.fcalc*1e9;

    % Postprocessing & do the plots
    freq = linspace(f_start,f_stop,201);
    port = calcPort(port, Sim.Sim_Path, freq);

    Zin = port.uf.tot ./ port.if.tot;
    s11 = port.uf.ref ./ port.uf.inc;

    % 1) Reflection coefficient (S11) plots
    % S11 absolute value
    S11_figure = figure('position',[600,100,900,900]);
    plot(freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2);
    xlim([Sim.fmin Sim.fmax]);
    ylim([-40 0]);
    set(gca, "linewidth",2, "fontsize", 14);
    grid on
    title('Reflection Coefficient S_{11}', 'FontSize', 16);
    xlabel('Frequency (GHz)','FontSize', 14);
    ylabel('Reflection Coefficient |S_{11}| (dB)','FontSize', 14);
    drawnow

    S11_output = strcat(Sim.output_path, sprintf('/svg_plots/%d_S11_output.svg', Sim.horn_number));
    print(S11_figure, S11_output, "-dsvg", "-Sxsize=900");

    S11_output = strcat(Sim.output_path, sprintf('/jpg_plots/%d_S11_output.jpg', Sim.horn_number));
    print(S11_figure, S11_output, "-djpg", "-Sxsize=900");

    % S11 Smith chart
    S11_smith_figure = figure('position',[600,100,900,900]);
    plotRefl(port, 'fmarkers', [f_start, f_stop], 'linewidth', 2, 'legend_fontsize', 20, 'markersize', 15);
    set(gca,"linewidth",2, "fontsize", 14, "XTick", -180:30:180, "YTick", -30:5:40);
    title('Reflection Coefficient S_{11}', 'FontSize', 16);
    drawnow

    S11_smith_output = strcat(Sim.output_path, sprintf('/svg_plots/%d_S11_smith_output.svg', Sim.horn_number));
    print(S11_smith_figure, S11_smith_output, "-dsvg", "-Sxsize=900");

    S11_smith_output = strcat(Sim.output_path, sprintf('/jpg_plots/%d_S11_smith_output.jpg', Sim.horn_number));
    print(S11_smith_figure, S11_smith_output, "-djpg", "-Sxsize=900");

    %%% NFFF plots

    % Calculate the far field at phi=0, 45 and at phi=90 degrees
    thetaRange = (0:0.2:359) - 180;
    disp('calculating far field at phi=[0 45 90] deg...');
    nf2ff = CalcNF2FF(nf2ff, Sim.Sim_Path, f0, thetaRange*pi/180, [0 45 90]*pi/180);

    Dlog=10*log10(nf2ff.Dmax);      % Calculate maximum Directivity in dB

    % Display some antenna parameters from above calculations
    disp(['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
    disp(['directivity: Dmax = ' num2str(Dlog) ' dBi']);

    % 2) Directivity
    farfield_directivity_figure = figure('position',[600,100,900,900]);
    plotFFdB(nf2ff, 'xaxis', 'theta', 'param', [1 2 3], 'linewidth', 2, 'legend_fontsize', 20);
    ylim([-30 25]);
    xlim([-180 180]);
    grid on
    set(gca,"linewidth",2, "fontsize", 14, "XTick", -180:30:180, "YTick", -30:5:40)
    title(sprintf('Farfield Directivity @ %.2f GHz',Sim.fcalc),'FontSize', 16);
    xlabel('Theta (degrees)','FontSize', 14);
    ylabel('Directivity (dBi)','FontSize', 14);
    drawnow

    farfield_directivity_output = strcat(Sim.output_path, sprintf('/svg_plots/%d_Farfield_directivity_output.svg', Sim.horn_number));
    print(farfield_directivity_figure, farfield_directivity_output, "-dsvg", "-Sxsize=900");

    farfield_directivity_output = strcat(Sim.output_path, sprintf('/jpg_plots/%d_Farfield_directivity_output.jpg', Sim.horn_number));
    print(farfield_directivity_figure, farfield_directivity_output, "-djpg", "-Sxsize=900");

    % 3) Plot Ludwig3 cross polar
    farfield_directivity_ludwig_figure = figure('position',[600,100,900,900]);
    plotFFcocx(nf2ff, 'xaxis', 'theta', 'param', [2], 'linewidth', 2, 'legend_fontsize', 20);
    ylim([-30 25]);
    xlim([-180 180]);
    grid on
    set(gca,"linewidth",2, "fontsize", 14, "XTick", -180:30:180, "YTick", -30:5:40)
    title(sprintf('Farfield Directivity with Ludwig3 XPOL @ %.2f GHz',Sim.fcalc),'FontSize', 16);
    xlabel('Theta (degrees)','FontSize', 14);
    ylabel('Directivity (dBi)','FontSize', 14);
    drawnow
    
    farfield_directivity_ludwig_output = strcat(Sim.output_path, sprintf('/svg_plots/%d_Farfield_directivity_ludwig_output.svg', 
                                                                                                    Sim.horn_number));
    print(farfield_directivity_ludwig_figure, farfield_directivity_ludwig_output, "-dsvg", "-Sxsize=900");

    farfield_directivity_ludwig_output = strcat(Sim.output_path, sprintf('/jpg_plots/%d_Farfield_directivity_ludwig_output.jpg', 
                                                                                                    Sim.horn_number));
    print(farfield_directivity_ludwig_figure, farfield_directivity_ludwig_output, "-djpg", "-Sxsize=900");

    % 4) Polar plot
    farfield_directivity_polar_figure = figure('position',[600,100,900,900]);
    leg=[];   %legend
    polarFF(nf2ff,'xaxis','theta','param',[1 2 3],'logscale',[-30 35], 'xtics', 12, 'linewidth', 2, 'legend_fontsize', 20);
    title(sprintf('Farfield Directivity @ %.2f GHz',Sim.fcalc),'FontSize', 16);
    xlabel('Theta (degrees)','FontSize', 14);
    ylabel('Directivity (dBi)','FontSize', 14);
    drawnow

    farfield_directivity_polar_output = strcat(Sim.output_path, sprintf('/svg_plots/%d_Farfield_directivity_polar_output.svg', 
                                                                                                    Sim.horn_number));
    print(farfield_directivity_polar_figure, farfield_directivity_polar_output, "-dsvg", "-Sxsize=900");

    farfield_directivity_polar_output = strcat(Sim.output_path, sprintf('/jpg_plots/%d_Farfield_directivity_polar_output.jpg', 
                                                                                                    Sim.horn_number));
    print(farfield_directivity_polar_figure, farfield_directivity_polar_output, "-djpg", "-Sxsize=900");

    %% Calculate 3D pattern
    %phiRange = sort(unique([-180:5:-100 -100:2.5:-50 -50:1:50 50:2.5:100 100:5:180]));
    %thetaRange = sort(unique([0:1:50 50:2:100 100:5:180]));
    phiRange = sort(unique([-180:1:-100 -100:1:-50 -50:1:50 50:1:100 100:1:180]));
    thetaRange = sort(unique([0:1:50 50:1:100 100:1:180]));

    disp('calculating 3D far field...');
    nf2ff = CalcNF2FF(nf2ff, Sim.Sim_Path, f0, thetaRange*pi/180, phiRange*pi/180, 'Verbose',2,'Outfile','nf2ff_3D.h5');

    radiation_pattern_figure = figure('position',[600,100,900,900]);
    colormap jet;
    plotFF3D(nf2ff, 'logscale', -40);        % plot 3D far field in dB

    radiation_patern_output = strcat(Sim.output_path, sprintf('/jpg_plots/%d_radiation_patern_output.jpg', Sim.horn_number));
    print(radiation_pattern_figure, radiation_patern_output, "-djpg", "-Sxsize=900");

    % Save far field in VTK to plot in ParaView
    E_far_normalized = nf2ff.E_norm{1}/max(nf2ff.E_norm{1}(:));
    DumpFF2VTK([Sim.Sim_Path '/Farfield.vtk'],E_far_normalized,thetaRange,phiRange,'scale', 0.008, 'logscale', -30, 
                                                                                                    'maxgain', Dlog);
endfunction