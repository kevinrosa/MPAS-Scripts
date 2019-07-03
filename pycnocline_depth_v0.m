% Kevin Rosa
% July 3, 2019
%
% similar to ssh_averaged_basin_v0.m, this script will average data on
% native mpas grid (will not interpolate to regular grid)

addpath(genpath('.'))

%%
type = 'mpaso.hist.am.timeSeriesStatsMonthly';

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 2:6;%:30;
run(i).color = rgb('red');
% i = i+1;
% run(i).name = 'High-resolution G case';
% run(i).short_name = 'high-res';
% run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
% run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
% run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
% run(i).years = 1:4;%:19;
% run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = '20180305.GM600.T62_oECv3.eos';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).years = 1; %25:30;%20:36;
run(i).color = rgb('blue');

%%
lon_range = [-130, -122];
lat_range = [38 42];

FIELDS3D = {'timeMonthly_avg_potentialDensity'};


for i = 1:length(run)
    files = {};
    for year = run(i).years
        dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.timeSeriesStatsMonthly.%04i-*-01.nc',year)));
        files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
    end

    [all.lon, all.lat] = read_mesh_file_lonlat(run(i).mesh_fi);
    all.areaCell = ncread(run(i).mesh_fi,'areaCell');
    all.H = ncread(run(i).mesh_fi,'bottomDepth');
    
    xy_inds = all.lon>lon_range(1) & all.lon<lon_range(end) & all.lat>lat_range(1) & all.lat<lat_range(end);
    
    % only keep cells in the lon lat range
    run(i).lon = all.lon(xy_inds);
    run(i).lat = all.lat(xy_inds);
    run(i).areaCell  = all.areaCell(xy_inds);
    run(i).H = all.H(xy_inds);
    
    z_level_depths = ncread(run(i).mesh_fi, 'refBottomDepth');
    depth_inds = find(z_level_depths < max(run(i).H(:)));
    depth_inds = [depth_inds; depth_inds(end)+1];

    run(i).z = z_level_depths(depth_inds);
    
    
    % Initialize matrices with NaNs:
    t_ind = 1;   % netcdf time index to read
    t_length = length(files);  % number of time indices 

    nans1d = NaN(t_length, 1);
    run(i).time = nans1d;
    
    nans3d = NaN([t_length, length(depth_inds), length(find(xy_inds))]);
    for F = FIELDS3D
        run(i).(F{1}) = nans3d;
    end

    tt = 1;
    for m = 1:length(files)
        data_fi = files{m};
        run(i).time(tt) = mpas_time(data_fi, t_ind, 'xtime_startMonthly');   
        
        % read data fields
        for F = FIELDS3D
            data = ncread(data_fi, F{1}, [depth_inds(1),1,t_ind], [length(depth_inds),Inf,1]);
            run(i).(F{1})(tt,:,:) = data(:, xy_inds);
            
        end
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
        tt = tt+1;
    end
    
    % turn zeros into NaNs
    for F = FIELDS3D
        run(i).(F{1})(run(i).(F{1})==0) = NaN;
    end

end

%%
i = 1;

rho_target = 1026;

t = 10;
figure
for j = 1:length(run(i).timeMonthly_avg_potentialDensity(t,1,:))
    line(run(i).timeMonthly_avg_potentialDensity(t,:,j), run(i).z)
end

set(gca,'ydir','reverse')

%%
run(i).pycnocline = NaN([length(run(i).time), length(run(i).lon)]);

for t = 1:length(run(i).time)
    for j = 1:length(run(i).timeMonthly_avg_potentialDensity(t,1,:))
        vert_inds = ~isnan(run(i).timeMonthly_avg_potentialDensity(t,:,j));
        run(i).pycnocline(t,j) = interp1(run(i).timeMonthly_avg_potentialDensity(t,vert_inds,j), run(i).z(vert_inds), rho_target);
        
    end
end

%% calculate month
dv = datevec(run(i).time);
run(i).month = dv(:,2);

%%
space_mean = mean(run(i).pycnocline, 2, 'omitnan');
month_mean = NaN(12,1);
for month = 1:12
    t_inds = run(i).month == month;
    run(i).month_mean(month) = mean(space_mean(t_inds));
    run(i).month_std(month) = std(space_mean(t_inds));
end

%%
figure
months = 1:12;
month_order = [8:12, 1:7];
time = [datenum(1,8:12,1), datenum(2,1:7,1)];
line(time, run(i).month_mean(month_order))
line(time, run(i).month_mean(month_order)+run(i).month_std(month_order), 'linestyle','--')
line(time, run(i).month_mean(month_order)-run(i).month_std(month_order), 'linestyle','--')


datetick('x','mmm')

set(gca,'ydir','reverse')
