% Kevin Rosa
% July 3, 2019
%
% similar to ssh_averaged_basin_v0.m, this script will average data on
% native mpas grid (will not interpolate to regular grid)
%
% TO DO:
% - Add area-weighted averaging for pycnocline depth
% - Only read xy_inds between the first and the last index

addpath(genpath('.'))

%%
type = 'mpaso.hist.am.timeSeriesStatsMonthly';

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 10:14;%:30;
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).years = 10:13;
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
% lon_range = [-130, -122];
% lat_range = [38 42];
lon_range = [-125, -123];
lat_range = [39 41];

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
    
    xy_inds = find(all.lon>lon_range(1) & all.lon<lon_range(end) & all.lat>lat_range(1) & all.lat<lat_range(end));
    
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
    
    nans3d = NaN([t_length, length(depth_inds), length(xy_inds)]);
    for F = FIELDS3D
        run(i).(F{1}) = nans3d;
    end

    tt = 1;
    for m = 1:length(files)
        data_fi = files{m};
        run(i).time(tt) = mpas_time(data_fi, t_ind, 'xtime_startMonthly');   
        
        % read data fields
        for F = FIELDS3D
            start = xy_inds(1);
            count = xy_inds(end) - start + 1;
            data = ncread(data_fi, F{1}, [depth_inds(1),start,t_ind], [length(depth_inds),count,1]);
            run(i).(F{1})(tt,:,:) = data(:, xy_inds-start+1);
            
        end
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
        tt = tt+1;
    end
    
    % turn zeros into NaNs
    for F = FIELDS3D
        run(i).(F{1})(run(i).(F{1})==0) = NaN;
    end

end


%% compare density profiles
target_lon = -126;
target_lat = 40;
t = 7;

figure
for i = 1:length(run)
    [~,ind] = min((run(i).lon-target_lon).^2 + ((run(i).lat-target_lat).*cosd(run(i).lat)).^2);
    
    line(run(i).timeMonthly_avg_potentialDensity(t,:,ind), run(i).z, 'color',run(i).color)
    
end

set(gca,'ydir','reverse')
ylim([-1 300])


%% new sub-region
% % run1 = run;
% lon_range = [-125, -123];
% lat_range = [39 41];
% 
% for i = 1:length(run)
%     xy_inds = find(run1(i).lon>lon_range(1) & run1(i).lon<lon_range(end) & run1(i).lat>lat_range(1) & run1(i).lat<lat_range(end));
%     
%     for F = FIELDS3D
%         start = xy_inds(1);
%         count = xy_inds(end) - start + 1;
%         data = ncread(data_fi, F{1}, [depth_inds(1),start,t_ind], [length(depth_inds),count,1]);
%         run(i).(F{1}) = run1(i).(F{1})(:,:,xy_inds);
% 
%     end
% end

%% solving a weird issue with the low-res run
% indices beyond depth range are a low density, not 0 or NaN
i = 3;
min_rho = 1010;
run(i).timeMonthly_avg_potentialDensity(run(i).timeMonthly_avg_potentialDensity<min_rho) = NaN;


%% calculate grid cell widths
i = 1;
lon_vec = -150:0.1:-110;
lat_vec = 20:60;
[LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
[~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, LON(:,1), LAT(1,:), 0);
run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;  % km

% generate transition line 
transition_contour = 15;
[T.lon, T.lat] = streamline_coords(LON, LAT, run(i).widthCell, transition_contour);

%% show region
figure
lon_range = [-130, -122];
lat_range = [38 42];

m_proj('lambert','long',[-140 -115],'lat',[30 50]);
% m_proj('lambert','long',[-125.5 -122.5],'lat',[38.5 41.5]);
hold on

m_gshhs_h('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')

m_line([lon_range(1),lon_range(1),lon_range(2),lon_range(2),lon_range(1)], [lat_range(1),lat_range(2),lat_range(2),lat_range(1),lat_range(1)], 'color','k','linewidth',3); 

m_line(T.lon, T.lat, 'color','k','linestyle','--','linewidth',2)

%%
set(gcf,'color','w')
saveas(gcf,'figures/upwelling/isopycnal_region_v01.png')

%%
i = 1;

rho_target = 1026;

t = 1;
figure
for j = 1:length(run(i).timeMonthly_avg_potentialDensity(t,1,:))
    line(run(i).timeMonthly_avg_potentialDensity(t,:,j), run(i).z)
end

set(gca,'ydir','reverse')

%% calculate pycnocline depth for each cell at each time
for i = 1:length(run)
    run(i).pycnocline = NaN([length(run(i).time), length(run(i).lon)]);

    for t = 1:length(run(i).time)
        for j = 1:length(run(i).timeMonthly_avg_potentialDensity(t,1,:))
            vert_inds = ~isnan(run(i).timeMonthly_avg_potentialDensity(t,:,j));
            run(i).pycnocline(t,j) = interp1(run(i).timeMonthly_avg_potentialDensity(t,vert_inds,j), run(i).z(vert_inds), rho_target);

        end
    end
end

%% calculate spatial mean pycnocline depth
for i = 1:length(run)
    dv = datevec(run(i).time);
    run(i).month = dv(:,2);

    space_mean = mean(run(i).pycnocline, 2, 'omitnan');
    
    percentiles = prctile(run(i).pycnocline, [25 50 75], 2);
    
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

% month_order = 1:12;
% time = datenum(1,1:12,1);

for i = 1:length(run)
    L(i) = line(time, run(i).month_mean(month_order), 'color',run(i).color,'linewidth',2);
    line(time, run(i).month_mean(month_order)+run(i).month_std(month_order), 'linestyle','--', 'color',run(i).color)
    line(time, run(i).month_mean(month_order)-run(i).month_std(month_order), 'linestyle','--', 'color',run(i).color)
    
%     line(time, run(i).percentiles(month_order, 2), 'color',run(i).color,'linewidth',2)

end
lg = legend(L, {run(:).short_name});



datetick('x','mmm')

set(gca,'ydir','reverse')
xlim([time(1) time(end)])

title(sprintf('%.1f isopycnal mean', rho_target))

%%
set(gcf,'color','w')
saveas(gcf,'figures/upwelling/isopycnal_months_smallerregion_wvariance_v02.png')

%%
i = 1;
t = 7;
figure
histogram(run(i).pycnocline(t,:))

%% plot spatial distribution of pycnocline depth — should show upwelling filaments for high res runs…
i = 3;
t = 7;

dx = 0.1;
lon_vec = lon_range(1):0.1:lon_range(2);
lat_vec = lat_range(1):0.1:lat_range(2);
[LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);

F = scatteredInterpolant(run(i).lon, run(i).lat, run(i).pycnocline(t,:)', 'linear','none'); 
FIELD = F(LON, LAT);

%%
figure
pcolor(LON, LAT, FIELD); shading flat

caxis([0 140])
