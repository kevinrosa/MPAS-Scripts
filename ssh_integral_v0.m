% Kevin Rosa
% June 3, 2019

%%

addpath(genpath('.'))

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).dir = '/scratch/kanga/runs/GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01/mpaso.hist.am.highFrequencyOutput/';
run(i).mesh_fi = '/scratch/kanga/runs/GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01/mpaso.rst.0001-08-01_00000.nc';
run(i).years = 2;%2:22;
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).dir = '/scratch/kanga/runs/20180208.GMPAS-IAF.T62_oRRS18v3.anvil/mpaso.hist.am.highFrequencyOutput/';
run(i).mesh_fi = '/scratch/kanga/runs/20180208.GMPAS-IAF.T62_oRRS18v3.anvil/ocean.oRRS18to6v3.scrip.181106.nc';
run(i).years = 20;%%20:36;
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).dir = '/scratch/kanga/runs/20180305.GM600.T62_oECv3.eos/mpaso.hist.am.highFrequencyOutput/';
run(i).mesh_fi = '/scratch/kanga/runs/20180305.GM600.T62_oECv3.eos/mpaso.rst.0050-01-01_00000.nc';
run(i).years = 60;%20:36;
run(i).color = rgb('blue');

%% download SSH

dx = 0.2;
lon_vec = -81:dx:-11;
lat_vec = 25:dx:50;

for i = 1:length(run)
    files = {};
    for year = run(i).years
        dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.highFrequencyOutput.%04i-*',year)));
        files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
    end

    t_ind = 1;
    t_length = length(files);  % number of time indices 

%     run(i).dist = DIST;
%     nans = NaN(t_length, length(DIST));
%     run(i).flux = nans;
%     run(i).lon = nans;
%     run(i).lat = nans;
%     run(i).time = NaN(t_length, 1);
    
%     run(i).speed = NaN([t_length, length(lon_vec), length(lat_vec)]);

    tt = 1;
    for m = 1%:length(files)

        [LON, LAT, ssh] = mpas_to_lonlat_meshgrid('ssh', run(i).mesh_fi, files{m}, lon_vec, lat_vec, t_ind);
%         [LON, LAT, ke_mat] = mpas_to_lonlat_meshgrid('kineticEnergyAtSurface', run(i).mesh_fi, files{m}, lon_vec, lat_vec, t_ind);
%         speed_mat = sqrt(2 * ke_mat);
%         run(i).speed(m,:,:) = speed_mat;

%         [C.lon, C.lat] = streamline_coords(LON, LAT, ssh, run(i).ssh_contour);
        
%         if ~isempty(C.lon)
%             [left, right, angle] = cross_stream_transects(C.lon, C.lat, width);
% 
%             flux = NaN(size(left.lon));
%             for k = 1:length(left.lon)
%                 flux(k) = integrate_along_transect(LON,LAT,speed_mat,left.lon(k),left.lat(k),right.lon(k),right.lat(k),N_integration);
%             end
% 
%             dist = distance_along_stream(C.lon, C.lat);
% 
%             run(i).flux(tt,:) = interp1(dist, flux, DIST);
%             run(i).angle(tt,:) = interp_angle(dist, angle, DIST);
%             run(i).lon(tt,:) = interp1(dist, C.lon, DIST);
%             run(i).lat(tt,:) = interp1(dist, C.lat, DIST);
%         end

        time_str = ncread(files{m},'xtime',[1,t_ind],[Inf,1]);
        run(i).time(tt) = datenum(time_str','yyyy-mm-dd');
        
        run(i).ssh = ssh;

        tt = tt+1;
    end
end

%% compute a land mask
i = 2;
[~,~,temp] = mpas_to_lonlat_meshgrid('temperatureAt250m', run(i).mesh_fi, fullfile(run(i).dir,'mpaso.hist.am.highFrequencyOutput.0020-01-01_00.00.00.nc'), lon_vec, lat_vec, t_ind);
mask = ones(size(temp));
mask(temp<-10) = NaN;

%%
figure
i = 2;
pcolor(LON, LAT, run(i).ssh.*mask); shading flat

%%
for i = 1:length(run)
    nansum(run(i).ssh(:).*mask(:))
end

%%
figure
set(gcf,'color','w')
row = 50;
for i = 1:length(run)
    line(LON(:,row), run(i).ssh(:,row),'color',run(i).color,'linewidth',2)
end
title(sprintf('SSH across %.1f ^o N',LAT(1,row)))
xlabel('Longitude (^o East)')
xlim([-80 -10])
legend(run(:).short_name)

set(gca,'fontsize',14)

% saveas(gcf, sprintf('figures/geoid/ssh_section01_v0.png'))

%%
nanmean(run(i).ssh(:,row) .* mask(:,row))
nanmean(run(i).ssh(:) .* mask(:))









