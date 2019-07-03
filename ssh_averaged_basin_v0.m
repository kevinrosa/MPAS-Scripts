addpath(genpath('.'))


i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = '20180305.GM600.T62_oECv3.eos';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).color = rgb('blue');

%% get ssh mean for all available files for each run so can smooth in time
% lon_range = [-81, -10];
% lat_range = [25 50];

lon_range = [-157, -115];
lat_range = [29 46];

for i = 1:length(run)
    dd = dir(fullfile(run(i).dir, 'mpaso.hist.am.highFrequencyOutput*'));
    files = fullfile({dd(:).folder}, {dd(:).name})';

    t_ind = 1;
    t_length = length(files);  % number of time indices 

    run(i).time = NaN(t_length, 1);
    run(i).ssh_mean = run(i).time; 

    tt = 1;
    for m = 1:length(files)
        
        data_fi = files{m};
        
        run(i).ssh_mean(tt) = mean_mpas_area_weighted('ssh', run(i).mesh_fi, data_fi, lon_range, lat_range, t_ind);
        
        run(i).time(tt) = mpas_time(data_fi, t_ind);
        
        tt = tt+1;
    end
end    

%% plot results 
cutoff_period_years = 5;
cutoff_period = cutoff_period_years * 364.25 * 24;  % in hours

for i = 1:length(run)
    figure
    set(gcf,'color','w')
    
    line(run(i).time, run(i).ssh_mean,'color',run(i).color)
    
    [run(i).ssh_mean_lp, run(i).time_lp] = lowpassfilter(run(i).time, run(i).ssh_mean, cutoff_period);
    
    line(run(i).time_lp, run(i).ssh_mean_lp, 'color','k','linewidth',2)
    
    datetick('x','yy')
    
%     set(gca,'fontsize',14)
    xlabel('Year')
    ylabel('Area-averaged SSH (m)')
    title(sprintf('%s average SSH with %i year lowpass', run(i).name, cutoff_period_years))
    
%     save_name = sprintf('figures/ssh_trend/ssh_atlantic_lon%.1fto%.1f_lat%.1fto%.1f_%s_v0.png', lon_range, lat_range,run(i).short_name);
%     saveas(gcf, save_name)
end

%% save data
for i = 1:length(run)
    time = run(i).time;
    ssh_mean = run(i).ssh_mean;
    
    save_name = sprintf('ssh_mean_highfreq_lon%.1fto%.1f_lat%.1fto%.1f_%s_v0.mat', lon_range, lat_range, run(i).code);
    save(save_name, 'time', 'ssh_mean')
end
for i = 1:length(run)
    time = run(i).time_lp;
    ssh_mean = run(i).ssh_mean_lp;
    
    save_name = sprintf('ssh_mean_lowpass%iyears_lon%.1fto%.1f_lat%.1fto%.1f_%s_v0.mat',cutoff_period_years, lon_range, lat_range, run(i).code);
    save(save_name, 'time', 'ssh_mean')
end

