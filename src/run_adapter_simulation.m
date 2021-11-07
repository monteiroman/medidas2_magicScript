    
function run_adapter_simulation(Sim, port, freq)
    %%
    %%
    %%  Parameters:
    %%          Sim: Horn simulation parameters.
    %%          port: Waveguide port for later simulation.
    %%
    %%  Returns:
    %%          void
        
    RunOpenEMS(Sim.Sim_Path, Sim.Sim_CSX, '--dump-statistics')

    port = calcPort( port, Sim.Sim_Path, freq);

    % must correct s21 by ratio of each port impedance
    % impedance of the coax is still not perfect due to mesh discretization 
    % in a cartesian grid

    s11 = port{1}.uf.ref./ port{1}.uf.inc;
    s21 = sqrt(real(port{1}.ZL_ref)./port{2}.ZL_ref).*port{2}.uf.ref./port{1}.uf.inc;

    %% plot s-parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S11_figure = figure('position',[600,100,900,900]);
    plot(freq*1e-9,20*log10(abs(s11)),'k-','Linewidth',2);
    xlim([freq(1) freq(end)]*1e-9);
    grid on;
    hold on;
    plot(freq*1e-9,20*log10(abs(s21)),'r--','Linewidth',2);
    set(gca, "linewidth",2, "fontsize", 14);
    l = legend('S_{11}','S_{21}','Location','southeast');
    set(l,'FontSize',12);
    title('S_{11} and S_{21}', 'FontSize', 16);
    ylabel('S-Parameter (dB)','FontSize',14);
    xlabel('frequency (GHz) \rightarrow','FontSize',14);

    S11_figure_output = strcat(Sim.output_path, sprintf('/jpg_plots/%d_S11_output.jpg', Sim.adapt_number));
    print(S11_figure, S11_figure_output, "-djpg", "-Sxsize=900");

    S11_figure_output = strcat(Sim.output_path, sprintf('/svg_plots/%d_S11_output.svg', Sim.adapt_number));
    print(S11_figure, S11_figure_output, "-dsvg", "-Sxsize=900");

    
    %% plot impedance of coax port
    Impedance_figure = figure('position',[600,100,900,900]);
    plot(freq*1e-9,real(port{1}.ZL_ref),'k-','Linewidth',2);
    xlim([freq(1) freq(end)]*1e-9);
    set(gca, "linewidth",2, "fontsize", 14);
    grid on;
    hold on;
    title('Impedance', 'FontSize', 16);
    set(l,'FontSize',12);
    ylabel('coax impedance','FontSize',14);
    xlabel('frequency (GHz) \rightarrow','FontSize',14);

    Impedance_figure_output = strcat(Sim.output_path, sprintf('/jpg_plots/%d_Impedance_output.jpg', Sim.adapt_number));
    print(Impedance_figure, Impedance_figure_output, "-djpg", "-Sxsize=900");

    Impedance_figure_output = strcat(Sim.output_path, sprintf('/svg_plots/%d_Impedance_output.svg', Sim.adapt_number));
    print(Impedance_figure, Impedance_figure_output, "-dsvg", "-Sxsize=900");

    %% Plot the field dumps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Great diagnostic tool
    % figure
    % dump_file = [Sim_Path '/Et.h5'];
    % PlotArgs.slice = {a/2*unit (b/2)*unit (.05)*unit};
    % PlotArgs.pauseTime=0.01;
    % PlotArgs.component=0;
    % PlotArgs.Limit = 'auto';
    % PlotHDF5FieldData(dump_file, PlotArgs);
endfunction