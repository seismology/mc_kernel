function plot_test_results(bool_plot, descriptions)

rundirs = dir('test_*');

rundirs = rundirs(bool_plot);

clock(1).desc = 'FFT routines';
clock(2).desc = 'Reading bwd field';
clock(3).desc = 'Reading fwd field';
clock(4).desc = 'KD-tree lookup (only fwd)';
clock(5).desc = 'Load_strain (fwd and bwd)';
clock(6).desc = 'NetCDF routines';
clock(7).desc = 'Rotate fields';
clock(8).desc = 'Buffer routines';
clock(9).desc = 'Monte Carlo routines';
clock(10).desc = 'Filtering and convolution';
clock(11).desc = 'Inversion mesh routines';
clock(12).desc = 'Kernel routines';
clock(13).desc = 'Initialization per task';
clock(14).desc = 'MPI communication with Master';

clock(1).short = 'fft';
clock(2).short = 'read_bwd';
clock(3).short = 'read_fwd';
clock(4).short = 'kdtree';
clock(5).short = 'load_strain';
clock(6).short = 'netcdf';
clock(7).short = 'rotate';
clock(8).short = 'buffer';
clock(9).short = 'montecarlo';
clock(10).short = 'filt_and_conv';
clock(11).short = 'inv_mesh';
clock(12).short = 'kernel';
clock(13).short = 'init_per_task';
clock(14).short = 'mpi_comm';

nruns = length(rundirs);
nclock = 14;
for irun = 1:nruns
    result = read_kerner_timing(rundirs(irun).name);
    
    for iclock =1:nclock
        clock(iclock).ncalls(:, irun)      = result.ncalls(:,iclock);
        clock(iclock).timepercall(:, irun) = result.timepercall(:,iclock);
        clock(iclock).timetotal(:, irun)   = result.timetotal(:,iclock);
        clock(iclock).timeratio(:, irun)   = result.timeratio(:,iclock);
    end
  
end

%% Plot ratios of most important parts
hfig = figure('Visible', 'off');
hold on;
plot(mean(clock(1).timeratio*100,1), 'k',   'LineWidth', 2)  % FFT
plot(mean(clock(6).timeratio*100,1), 'b--', 'LineWidth', 2)  % NetCDF
plot(mean(clock(8).timeratio*100,1), 'b-.',   'LineWidth', 2)  % Buffer
plot(mean(clock(8).timeratio*100,1)+ ...
     mean(clock(6).timeratio*100,1), 'b',   'LineWidth', 2)  % Buffer
plot(mean(clock(10).timeratio*100,1), 'r',  'LineWidth', 2) % Filter and Convolution

legend({'FFT calls', 'NetCDF calls', 'Buffer calls', 'NetCDF + Buffer', 'Filtering and convolution'}, ...
        'Location', 'NorthWest')
set(gca, 'XTick', 1:nruns, 'XTickLabel', descriptions)
xlabel('Tests')
ylim([0, 100])
ylabel('ratio of total time in %')
title('Ratio of most important program parts')
fnam = sprintf('test_results_comparison_ratio');
print('-dpng', fnam)
close(hfig)

%% Total time of most important parts
hfig = figure('Visible', 'off');
hold on;
plot(sum(clock(1).timetotal,1), 'k',   'LineWidth', 2)  % FFT
plot(sum(clock(6).timetotal,1), 'b--', 'LineWidth', 2)  % NetCDF
plot(sum(clock(8).timetotal,1), 'b-.',   'LineWidth', 2)  % Buffer
plot(sum(clock(8).timetotal,1)+ ...
     sum(clock(6).timetotal,1), 'b',   'LineWidth', 2)  % Buffer
plot(sum(clock(10).timetotal,1), 'r',  'LineWidth', 2) % Filter and Convolution

legend({'FFT calls', 'NetCDF calls', 'Buffer calls', 'NetCDF + Buffer', 'Filtering and convolution'}, ...
        'Location', 'NorthWest')
set(gca, 'XTick', 1:nruns, 'XTickLabel', descriptions)
xlabel('Tests')
% ylim([0, 100])
ylabel('total time in s')
title('Total time of most important program parts')
fnam = sprintf('test_results_comparison_totaltime');
print('-dpng', fnam)
close(hfig)

%% More plots

for iclock = 1:nclock
    %% Time per call
    hfig = figure('Visible', 'off');
    hold on;
    plot(clock(iclock).timepercall', 'o')
    plot(mean(clock(iclock).timepercall,1), 'Linewidth', 2)
    set(gca, 'XTick', 1:nruns, 'XTickLabel', descriptions)
    xlabel('Tests')
    ylim([0, max(max(clock(iclock).timepercall))*1.1])
    ylabel('time per call / s')
    
    title(clock(iclock).desc)
    fnam = sprintf('test_results_%s_timepercall', clock(iclock).short);
    print('-dpng', fnam)
    close(hfig)
    
    %% Number of calls
    hfig = figure('Visible', 'off');
    hold on;
    plot(clock(iclock).ncalls', 'o')
    plot(mean(clock(iclock).ncalls,1), 'Linewidth', 2)
    set(gca, 'XTick', 1:nruns, 'XTickLabel', descriptions)
    xlabel('Tests')
    ylim([0, max(max(clock(iclock).ncalls))*1.1])
    ylabel('Number of calls')
    
    title(clock(iclock).desc)
    fnam = sprintf('test_results_%s_ncalls', clock(iclock).short);
    print('-dpng', fnam)
    close(hfig)
    
    %% Total time
    hfig = figure('Visible', 'off');
    hold on;
    plot(clock(iclock).timetotal', 'o')
    plot(mean(clock(iclock).timetotal,1), 'Linewidth', 2)
    set(gca, 'XTick', 1:nruns, 'XTickLabel', descriptions)
    xlabel('Tests')
    ylim([0, max(max(clock(iclock).timetotal))*1.1])
    ylabel('total time / s')
    
    title(clock(iclock).desc)
    fnam = sprintf('test_results_%s_timetotal', clock(iclock).short);
    print('-dpng', fnam)
    close(hfig)
    
    
    
    %% Time ratio
    hfig = figure('Visible', 'off');
    hold on;
    plot(clock(iclock).timeratio' * 100, 'o')
    plot(mean(clock(iclock).timeratio * 100, 1), 'Linewidth', 2)
    set(gca, 'XTick', 1:nruns, 'XTickLabel', descriptions)
    xlabel('Tests')
    ylim([0, 100])
    ylabel('ratio of total time in %')
    
    title(clock(iclock).desc)
    fnam = sprintf('test_results_%s_timeratio', clock(iclock).short);
    print('-dpng', fnam)
    close(hfig)
    
   
end



