% Kevin Rosa
% June 28, 2019

addpath(genpath('.'))

%%
load('ubar_structure_20190628.122449.mat')

%% time-means
for i = 1:length(run)    
    run(i).ubar_mean = squeeze(mean(run(i).timeMonthly_avg_velocityZonalDepthIntegrated,1));    
    run(i).vbar_mean = squeeze(mean(run(i).timeMonthly_avg_velocityMeridionalDepthIntegrated,1));
end


%%

type = 'mpaso.hist.am.timeSeriesStatsMonthly';

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 20;%20:30;
run(i).color = rgb('red');
run(i).levelp1 = 3;
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).years = 10;%10:19;
run(i).color = rgb('black');
run(i).levelp1 = 10;
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = '20180305.GM600.T62_oECv3.eos';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).years = 1; %25:30;%20:36;
run(i).color = rgb('blue');

%%

xrange = [-84 -40];
yrange = [0 50];
dx = 0.1;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

FIELDS = {'timeMonthly_avg_velocityZonalDepthIntegrated','timeMonthly_avg_velocityMeridionalDepthIntegrated','timeMonthly_avg_waterColumnThickness'};

for i = 1:length(run)
    files = {};
    for year = run(i).years
        dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.timeSeriesStatsMonthly.%04i-*ubar.nc',year)));
        files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
    end

    t_ind = 1;
    t_length = length(files);  % number of time indices 

    nans1d = NaN(t_length, 1);
    run(i).time = nans1d;
    
    nans3d = NaN([t_length, length(lon_vec), length(lat_vec)]);
    for F = FIELDS
        run(i).(F{1}) = nans3d;
    end

    tt = 1;
    for m = 1%:length(files)
        
        data_fi = files{m};
        
%         run(i).time(tt) = mpas_time(data_fi, t_ind, 'xtime_startMonthly');            

        
        % Calculate LON LAT and land-sea mask (only do this once)
        if tt == 1
            [LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
            run(i).LON = LON;
            run(i).LAT = LAT;
            
            run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);
            
        end
        
        % read data fields
        for F = FIELDS
            [~,~,data] = mpas_to_lonlat_meshgrid(F{1}, run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
            run(i).(F{1})(tt,:,:) = data .* run(i).mask;
        end
        
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
        tt = tt+1;
    end
end









%%
t = 1;
for i = 1:length(run)    
    run(i).ubar_mean = squeeze(run(i).timeMonthly_avg_velocityZonalDepthIntegrated(t,:,:));    
    run(i).vbar_mean = squeeze(run(i).timeMonthly_avg_velocityMeridionalDepthIntegrated(t,:,:));
end

%% calculate grid cell widths
for i = 1:length(run)
    [~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
    run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;  % km
end

%%
figure
i = 1;
transition_contour = 15;
contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)

%% generate transition line 
i = 1;
[T.lon, T.lat] = streamline_coords(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour);


%% view transition line
figure
xrange = [min(run(i).LON(:)), -40];%max(run(i).LON(:))];
yrange = [min(run(i).LAT(:)), max(run(i).LAT(:))];

m_proj('lambert','long',xrange,'lat',yrange);
hold on
m_line(T.lon, T.lat, 'color','k','linewidth',3); 

m_gshhs_l('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')

%% compute boundary-normal unit vectors
width = 100;
[left,right,angle] = cross_stream_transects(T.lon, T.lat, width);

[dx,dy] = lonlat_to_dxdy(T.lon(:), T.lat(:), right.lon(:), right.lat(:));

% east and north components of unit vector pointing out of coastal region 
T.ii = dx ./ width;
T.jj = dy ./ width;


%% interpolate and rotate velocities onto each boundary point
for i = 1:length(run)
    
    run(i).T = T;

    run(i).T.ubar_east = interp2(run(i).LON', run(i).LAT', run(i).ubar_mean', T.lon, T.lat);
    run(i).T.vbar_east = interp2(run(i).LON', run(i).LAT', run(i).vbar_mean', T.lon, T.lat);
    
    depth_mean = squeeze(mean(run(i).timeMonthly_avg_waterColumnThickness, 1,'omitnan'));
    run(i).T.depth = interp2(run(i).LON', run(i).LAT', depth_mean', T.lon, T.lat);

    % rotate to boundary-normal velocities
    run(i).T.ubar_normal = run(i).T.ubar_east(:) .* T.ii + run(i).T.vbar_east(:) .* T.jj;

end

%%
figure

for i = 1:length(run)
    line(T.lat, run(i).T.ubar_normal,'color',run(i).color)
end



%% spacing between points along boundary line
lat_ind_range = [24.5, 29];
lat_ind_range = [10, 15];
% lat_ind_range = [10, 45];
lat_ind_range = [15, 25]; 
lat_ind_range = [10, 30];
lat_ind_range = [0, 25];

inds = find( T.lat>lat_ind_range(1) & T.lat<lat_ind_range(2) );

% inds_w = [inds(1)-1, inds, inds(end)+1];  % indices used for width calculation 
% % centered-difference method
% dlon = (T.lon(inds_w(3:end)) - T.lon(inds_w(1:end-2))) ./ 2;
% dlat = (T.lat(inds_w(3:end)) - T.lat(inds_w(1:end-2))) ./ 2;
inds_w = [inds, inds(end)+1];  % indices used for width calculation 
% forward-difference method
dlon = (T.lon(inds_w(2:end)) - T.lon(inds_w(1:end-1)));
dlat = (T.lat(inds_w(2:end)) - T.lat(inds_w(1:end-1)));

width = (pi/180)*6371*1e3 * sqrt( (dlon.*cosd(T.lat(inds))).^2 + dlat.^2 );  % meters

% Transports
figure

for i = 1:length(run)

    transport = run(i).T.ubar_normal(inds)' .* run(i).T.depth(inds) .* width;
%     transport = run(i).T.ubar_east(inds) .* run(i).T.depth(inds) .* width;

    line(T.lat(inds), transport,'color',run(i).color)
    
    fprintf('%s: %f\n', run(i).short_name, sum(transport,'omitnan') * 1e-6)
    
end
% xlim([23, 32])

%%
i = 2;
figure
speed = sqrt(run(i).ubar_mean.^2 + run(i).vbar_mean.^2);

pcolor(run(i).LON, run(i).LAT, speed); shading flat
hold on
uu = run(i).ubar_mean./speed;
vv = run(i).vbar_mean./speed;
s = 4;
quiver(run(i).LON(1:s:end,1:s:end), run(i).LAT(1:s:end,1:s:end), uu(1:s:end,1:s:end), vv(1:s:end,1:s:end), 'color','w')

line(T.lon,T.lat,'color','r','linewidth',3)

%%
figure
contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)




%%
for k = 1:length(T.lon)
    line([T.lon(k) right.lon(k)], [T.lat(k) right.lat(k)], 'color','r')
end


%% T R A S H
%% at each point along contour, show normal vel
i = 2;

scale = 500 / 0.04;  % 0.04 m/s is shown as a 100 km displacement

displacement = run(i).T.ubar_normal * scale;  % km
dx = displacement .* T.ii;
dy = displacement .* T.jj;
[lon, lat] = dxdy_to_lonlat(dx(:), dy(:), T.lon(:), T.lat(:));

figure
line(T.lon, T.lat, 'color',0.4*[1,1,1],'linewidth',1)

line(lon, lat, 'color',run(i).color,'linewidth',2)

%%
figure
m_proj('lambert','long',xrange,'lat',yrange);
hold on
m_line(T.lon, T.lat, 'color',0.4*[1,1,1],'linewidth',2); 

m_gshhs_l('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')

i = 1;
m_line(lon, lat, 'color',run(i).color,'linewidth',3)

i = 2;
m_line(lon, lat, 'color',run(i).color,'linewidth',3)

