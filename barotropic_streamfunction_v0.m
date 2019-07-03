% Kevin Rosa
% June 28, 2019

%%
addpath(genpath('.'))

%%
type = 'mpaso.hist.am.timeSeriesStatsMonthly';

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 20:30;
run(i).color = rgb('red');
run(i).levelp1 = 3;
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).years = 10:19;
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

xrange = [-97 0];
yrange = [10 50];
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
    for m = 1:length(files)
        
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
save(sprintf('ubar_structure_%s.mat',datestr(now,'yyyymmdd.HHMMSS')), 'run','-v7.3')

%%

%% time-means
for i = 1:length(run)    
    run(i).ubar_mean = squeeze(mean(run(i).timeMonthly_avg_velocityZonalDepthIntegrated,1));    
    run(i).vbar_mean = squeeze(mean(run(i).timeMonthly_avg_velocityMeridionalDepthIntegrated,1));
end


%%
crange = [-1 1] * 0.1;
dc = 0.01;
bins = 2*crange(1):dc:crange(2)*4;
cmap = cmocean('balance',length(bins)-1);

for i = 1:2
    figure
    contourf(run(i).LON, run(i).LAT, run(i).vbar_mean, bins)
    
    colormap(cmap)
    caxis(crange)
    colorbar
end

%% transport through a section calculation
target_lat = 30;
lon_range = [-65 -20];

run_inds = [3,1,2];
for i = run_inds
    [~,eta] = min(abs(run(i).LAT(1,:)-target_lat));
    xi = run(i).LON(:,eta)>=lon_range(1) & run(i).LON(:,eta)<=lon_range(2);
    
    lon = run(i).LON(xi,eta);
    lat = run(i).LAT(xi,eta);

    [dx,~] = lonlat_to_dxdy(lon(1),lat(1),lon(2),lat(2));  % in km
    
    depth = squeeze(mean(run(i).timeMonthly_avg_waterColumnThickness, 1));
    
    run(i).trans = sum(run(i).vbar_mean(xi,eta) .* depth(xi,eta) * dx*1000) * 1e-6;  % Sverdrups
    
%     fprintf('%s mean: %.3f m/s \n', run(i).short_name, mean(run(i).v_mean(xi,eta)))    
    fprintf('%s transport: %.2f Sv \n', run(i).short_name, run(i).trans)    
    
end


%% grid spacing matrices DX and DY
LON = run(i).LON;
LAT = run(i).LAT;
[dx_vec,~] = lonlat_to_dxdy(LON(1,:), LAT(1,:), LON(2,:), LAT(2,:));
DX = repmat(dx_vec(:)', [length(LON(:,1)), 1]) * 1000;
[~,dy_vec] = lonlat_to_dxdy(LON(:,1), LAT(:,1), LON(:,2), LAT(:,2));
DY = repmat(dy_vec(:), [1,length(LON(1,:))]) * 1000;

%%
for i = 1:length(run)
    run(i).psi = NaN(size(run(i).ubar_mean));
    depth = squeeze(mean(run(i).timeMonthly_avg_waterColumnThickness, 1));

    for eta = 1:length(run(i).LON(1,:))

        run(i).psi(:,eta) = cumsum(run(i).vbar_mean(:,eta) .* depth(:,eta) .* DX(:,eta), 'omitnan') * 1e-6;

    end
end
%%
for i = 1:length(run)
figure(100+i)
clf
set(gcf,'position',[252 313 1021 562], 'color','w')

crange = [-10 60];
dc = 1;
% bins = 2*crange(1):dc:crange(2);
% bins = -80:dc:40;
bins = crange(1):dc:crange(2);
bins = -60:dc:60;

cmap =flipud(cbrewer('div','Spectral',length(bins)-1,'pchip'));

cont_lines = [0, 10, 20];

contourf(run(i).LON, run(i).LAT, run(i).psi .* run(i).mask, bins, 'linestyle','none'); 
hold on
contour(run(i).LON, run(i).LAT, run(i).psi .* run(i).mask, cont_lines, 'color','k');
colorbar

caxis(crange)
colormap(cmap)

save_name = sprintf('figures/streamfunction/bar_streamfun-v_%s_v01.png',run(i).code);
saveas(gcf,save_name)
end


%%













%%
i = 2;
figure

crange = [0 1] * 0.05;
dc = 0.001;
bins = 2*crange(1):dc:crange(2)*4;

contourf(run(i).LON, run(i).LAT, sqrt(run(i).ubar_mean.^2 + run(i).vbar_mean.^2), bins)

caxis(crange)



%%
i = 1;

t = 1;
dz = run(i).timeMonthly_avg_waterColumnThickness(t,:,:);
ubar = run(i).timeMonthly_avg_velocityZonalDepthIntegrated(t,:,:);
vbar = run(i).timeMonthly_avg_velocityMeridionalDepthIntegrated(t,:,:);

run(i).transport = dz .* sqrt(ubar.^2 + vbar.^2) * 1e-6;  % Sv

%%
figure
% contourf(run(i).LON, run(i).LAT, squeeze(run(i).transport(t,:,:)))
contourf(run(i).LON, run(i).LAT, squeeze(vbar))

colorbar

