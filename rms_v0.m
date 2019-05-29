% RMS of surface velocities
% 
% Kevin Rosa
% May 28, 2019

% [LON, LAT, FIELD] = mpas_to_lonlat_meshgrid(field_to_read, mesh_fi, data_fi, lon_vec, lat_vec, t_ind)

i = 1;
run(i).name = 'NA8';
run(i).dir = fullfile('/scratch/kanga/runs/standalone_hoch_mesh',run(i).name);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).color = 'k';
i = i+1;
run(i).name = 'CUSP8';
run(i).dir = fullfile('/scratch/kanga/runs/standalone_hoch_mesh',run(i).name);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).color = 'r';
i = i+1;
run(i).name = 'CUSP12';
run(i).dir = fullfile('/scratch/kanga/runs/standalone_hoch_mesh',run(i).name);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).color = rgb('orange');
i = i+1;
run(i).name = 'CUSP20';
run(i).dir = fullfile('/scratch/kanga/runs/standalone_hoch_mesh',run(i).name);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).color = rgb('green');
i = i+1;
run(i).name = 'CUSP30';
run(i).dir = fullfile('/scratch/kanga/runs/standalone_hoch_mesh',run(i).name);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).color = rgb('blue');

%% transects
j = 1;
transect(j).name = 'New Jersey';
transect(j).lon = [-69.2, -74.0];
transect(j).lat = [36.6, 40.6];
j = j+1;
transect(j).name = 'Cape Hatteras';
transect(j).lon = [-72.5, -75.3];
transect(j).lat = [35.2, 38.2];

%%
for i = 1:length(run)
    files = dir(fullfile(run(i).dir, '*high*.nc'));
    
    dx = 0.1;
    lon_vec = -84:dx:-54;
    lat_vec = 30:dx:42;
    sz = [length(lon_vec), length(lat_vec)];
    
    ke_mat = NaN([sz, length(files)]);  % initialize empty matrix
    
    t_ind = 1;
    
    for month = 1:length(files)
        data_fi = fullfile(files(month).folder, files(month).name);
        [LON, LAT, ke_mat(:,:,month)] = mpas_to_lonlat_meshgrid('kineticEnergyAtSurface', run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
    end
    speed_mat = sqrt(2 * ke_mat);
    
    run(i).lon = LON;
    run(i).lat = LAT;
    run(i).speed_rms = rms(speed_mat, 3);
end


%%
addpath(genpath('.'))

%%
save_dir = 'figures/rms_v0';

xrange = [-85 -55];
yrange = [lat_vec(1), lat_vec(end)];

FS = 14;
crange = [0 1.5];

rms_cont = 0.4;

%%
for i = 1:length(run)

figure(110+i)
clf
set(gcf,'color','w','position',[50 354 888 585])
m_proj('lambert','long',xrange,'lat',yrange);

m_pcolor(run(i).lon, run(i).lat, run(i).speed_rms)

hold on
m_contour(run(i).lon, run(i).lat, run(i).speed_rms, rms_cont*[1,1], 'color','w','linewidth',2)

m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
% m_grid('linestyle','none','linewidth',2,'tickdir','out','fontsize',14)
m_grid('linestyle','none','box','fancy','tickdir','out','fontsize',FS)

for j = 1:length(transect)
    m_line(transect(j).lon, transect(j).lat, 'color',rgb('red'),'linewidth',2,'marker','o','markerfacecolor','w')
end

caxis(crange)

title(run(i).name)

cb = colorbar;
ylabel(cb, 'Surface Speed RMS (m/s)')

set(gca,'fontsize',FS)

colormap('jet')

saveas(gcf, fullfile(save_dir, sprintf('speed_rms_pcolor_%s_showtransects_v0.png', run(i).name)))
end

%%
figure(210)
clf
set(gcf,'color','w','position',[50 354 888 585])
m_proj('lambert','long',xrange,'lat',yrange);

m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linestyle','none','box','fancy','tickdir','out','fontsize',FS)

hold on
for i = 1:2
    
    [~,handle(i)] = m_contour(run(i).lon, run(i).lat, run(i).speed_rms, rms_cont*[1,1], 'color',run(i).color,'linewidth',2);
    
end
title(sprintf('%.1f m/s RMS Contours: NA8 vs. CUSP8',rms_cont))

set(gca,'fontsize',FS)

% m_legend(handle, run(:).name)

saveas(gcf, fullfile(save_dir, sprintf('speed_rms_contours_2runs_v0.png')))

