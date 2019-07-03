% Kevin Rosa
% July 2, 2019

addpath(genpath('.'))

%%
load('ubar_structure_20190628.122449.mat')

%% time-means
for i = 1:length(run)    
    run(i).ubar_mean = squeeze(mean(run(i).timeMonthly_avg_velocityZonalDepthIntegrated,1));    
    run(i).vbar_mean = squeeze(mean(run(i).timeMonthly_avg_velocityMeridionalDepthIntegrated,1));
    run(i).depth_mean = squeeze(mean(run(i).timeMonthly_avg_waterColumnThickness,1));
end

%% calculate grid cell widths
for i = 1:length(run)
    [~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
    run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;  % km
end

%%
i = 2;
figure
clf
pcolor(run(i).LON, run(i).LAT, run(i).vbar_mean .* run(i).depth_mean); shading flat

crange = [-1 1] * 100;
caxis(crange)
cmap = flipud(cbrewer('div','RdBu',21,'pchip'));
colormap(cmap)

%% meridional transport through section
target_lat = 30;
% lon_range = [-82 -8];
lon_range = [-65 -10];

run_inds = [3,1,2];

for i = run_inds
    [~,eta] = min(abs(run(i).LAT(1,:)-target_lat));
    xi = run(i).LON(:,eta)>=lon_range(1) & run(i).LON(:,eta)<=lon_range(2);
    
    lon = run(i).LON(xi,eta);
    lat = run(i).LAT(xi,eta);

    [dx,~] = lonlat_to_dxdy(lon(1),lat(1),lon(2),lat(2));  % in km
    
%     run(i).trans = sum(run(i).vbar_mean(xi,eta) .* run(i).depth_mean(xi,eta) * dx*1000, 'omitnan') * 1e-6;  % Sverdrups
    run(i).trans = squeeze(run(i).vbar_mean(xi,eta) .* run(i).depth_mean(xi,eta) * dx*1000 * 1e-6);  % Sverdrups
    
    fprintf('%s transport: %.2f Sv \n', run(i).short_name, sum(run(i).trans,'omitnan'))    
    
end

%%
figure
for i = 1:length(run)
    line(lon, run(i).trans, 'color',run(i).color,'linewidth',2)
    
end

%% meridional transport through section **w/ std
target_lat = 30;
% lon_range = [-82 -8];
lon_range = [-65 -10];

run_inds = [3,1,2];

for i = run_inds
    [~,eta] = min(abs(run(i).LAT(1,:)-target_lat));
    xi = run(i).LON(:,eta)>=lon_range(1) & run(i).LON(:,eta)<=lon_range(2);
    
    lon = run(i).LON(xi,eta);
    lat = run(i).LAT(xi,eta);

    [dx,~] = lonlat_to_dxdy(lon(1),lat(1),lon(2),lat(2));  % in km
    
%     run(i).trans = sum(run(i).vbar_mean(xi,eta) .* run(i).depth_mean(xi,eta) * dx*1000, 'omitnan') * 1e-6;  % Sverdrups
    vel = run(i).timeMonthly_avg_velocityMeridionalDepthIntegrated(:,xi,eta);
    depth = run(i).timeMonthly_avg_waterColumnThickness(:,xi,eta);
    
    trans = sum(vel .* depth .* dx*1e3 * 1e-6, 2, 'omitnan');  % time-varying transport (Sv) 
            
    fprintf('%s transport: %.2f +/- %.3f Sv \n', run(i).short_name, mean(trans), std(trans))    
    
end








