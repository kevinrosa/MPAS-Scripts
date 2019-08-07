% Kevin Rosa
% July 8, 2019
%
% Will try the monthly climatology calculation again but use mixed layer
% depth from the the highfrequency output instead.
%  
% similar to ssh_averaged_basin_v0.m, this script will average data on
% native mpas grid (will not interpolate to regular grid)
%
% TO DO:
% - Add area-weighted averaging for pycnocline depth
% - Only read xy_inds between the first and the last index

addpath(genpath('.'))

%%
type = 'mpaso.hist.am.highFrequencyOutput';

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 10:19;%:30;
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).years = 30:39;
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = 'EC60to30_G_case';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).years = 50:60;
run(i).color = rgb('blue');

%%
lon_range = [-130, -122];
lat_range = [38 42];

FIELDS2D = {'dThreshMLD','tThreshMLD'};


for i = 1:length(run)
    files = {};
    for year = run(i).years
        dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.highFrequencyOutput.%04i-*',year)));
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
    
    
    % Initialize matrices with NaNs:
    t_ind = 1;   % netcdf time index to read
    t_length = length(files);  % number of time indices 

    nans1d = NaN(t_length, 1);
    run(i).time = nans1d;
    
    nans2d = NaN([t_length, length(find(xy_inds))]);
    for F = FIELDS2D
        run(i).(F{1}) = nans2d;
    end

    tt = 1;
    for m = 1:length(files)
        data_fi = files{m};
        run(i).time(tt) = mpas_time(data_fi, t_ind);   
        
        % read data fields
        for F = FIELDS2D
            data = ncread(data_fi, F{1}, [1,t_ind], [Inf,1]);
            run(i).(F{1})(tt,:) = data(xy_inds);
            
        end
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
        tt = tt+1;
    end
    
    % turn zeros into NaNs
    for F = FIELDS2D
        run(i).(F{1})(run(i).(F{1})==0) = NaN;
    end

end




%% calculate spatial mean MLD
F = 'dThreshMLD';

for i = 1:length(run)
    dv = datevec(run(i).time);
    run(i).month = dv(:,2);

    space_mean = mean(run(i).(F), 2, 'omitnan');
    
    percentiles = prctile(run(i).(F), [10 25 75], 2);
    
    month_mean = NaN(12,1);
    for month = 1:12
        t_inds = run(i).month == month;
        run(i).month_mean(month) = mean(space_mean(t_inds));
        run(i).month_std(month) = std(space_mean(t_inds));
        run(i).percentiles(month,:) = mean(percentiles(t_inds,:),1);
    end
end

%%
figure
months = 1:12;
month_order = [8:12, 1:7];
time = [datenum(1,8:12,1), datenum(2,1:7,1)];

month_order = 1:12;
time = datenum(1,1:12,1);

for i = 1:length(run)
%     line(time, run(i).month_mean(month_order), 'color',run(i).color,'linewidth',2)
%     line(time, run(i).month_mean(month_order)+run(i).month_std(month_order), 'linestyle','--', 'color',run(i).color)
%     line(time, run(i).month_mean(month_order)-run(i).month_std(month_order), 'linestyle','--', 'color',run(i).color)
%     
    line(time, run(i).percentiles(month_order, 1), 'color',run(i).color,'linewidth',2)

end
% legend(run(:).short_name)

datetick('x','mmm')

set(gca,'ydir','reverse')
xlim([time(1) time(end)])

% title(sprintf('%.1f isopycnal 25th percentile', rho_target))

%%
set(gcf,'color','w')
% saveas(gcf,'figures/upwelling/isopycnal_months_25th_v01.png')

%%
i = 1;
t = 7;
figure
histogram(run(i).tThreshMLD(t,:))

%% plot spatial distribution of pycnocline depth — should show upwelling filaments for high res runs…
i = 1;
t = 7;

dx = 0.1;
lon_vec = lon_range(1):0.1:lon_range(2);
lat_vec = lat_range(1):0.1:lat_range(2);
[LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);

F = scatteredInterpolant(run(i).lon, run(i).lat, run(i).tThreshMLD(t,:)', 'linear','none'); 
FIELD = F(LON, LAT);

%%
figure
pcolor(LON, LAT, FIELD); shading flat

caxis([0 40])
