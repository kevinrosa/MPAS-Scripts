% RMS of surface velocities
% 
% Kevin Rosa
% May 28, 2019

i = 1;
run(i).name = 'NA8';
run(i).dir = fullfile('/scratch/kanga/runs/standalone_hoch_mesh',run(i).name);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).color = rgb('blue');
i = i+1;
run(i).name = 'CUSP8';
run(i).dir = fullfile('/scratch/kanga/runs/standalone_hoch_mesh',run(i).name);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).color = rgb('orange');
i = i+1;
run(i).name = 'CUSP12';
run(i).dir = fullfile('/scratch/kanga/runs/standalone_hoch_mesh',run(i).name);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
% run(i).color = rgb('orange');
i = i+1;
run(i).name = 'CUSP20';
run(i).dir = fullfile('/scratch/kanga/runs/standalone_hoch_mesh',run(i).name);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
% run(i).color = rgb('green');
i = i+1;
run(i).name = 'CUSP30';
run(i).dir = fullfile('/scratch/kanga/runs/standalone_hoch_mesh',run(i).name);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
% run(i).color = rgb('blue');

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
save_dir = 'figures/rms_v1';

xrange = [-85 -55];
yrange = [lat_vec(1), lat_vec(end)];

FS = 14;
crange = [0 1.5];

rms_cont = 0.4;
white_cont = 0.4;

%% convoluted way to make colormap 
conts = crange(1):0.05:crange(end);
cmap = jet(length(conts)-1);
vals = (conts(1:end-1)+conts(2:end))/2;
inds = vals<=white_cont;
cmap(inds,:) = repmat([1,1,1],length(find(inds)),1);
cmap(~inds,:) = jet(length(find(~inds)));

%%
for i = 1:length(run)

figure(110+i)
clf
set(gcf,'color','w','position',[50 354 888 585])
m_proj('lambert','long',xrange,'lat',yrange);

m_pcolor(run(i).lon, run(i).lat, run(i).speed_rms)

hold on
% m_contour(run(i).lon, run(i).lat, run(i).speed_rms, rms_cont*[1,1], 'color','w','linewidth',2)

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

colormap(cmap)

saveas(gcf, fullfile(save_dir, sprintf('speed_rms_pcolor_%s_showtransects_v1.png', run(i).name)))
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

for j = 1:length(transect)
    m_line(transect(j).lon, transect(j).lat, 'color',rgb('red'),'linewidth',2,'marker','o','markerfacecolor','w')
end

title(sprintf('%.1f m/s RMS Contours: NA8 vs. CUSP8',rms_cont))

set(gca,'fontsize',FS)

% m_legend(handle, run(:).name)

saveas(gcf, fullfile(save_dir, sprintf('speed_rms_contours_2runs_showtransects_v1.png')))

%%
%
%
%
%
%% Gulf stream pathlines
% run(i).S will hold streamline info

years = 7:10;
ssh_contour = -0.2;
runs = 1:2;  % run indices to run calculation on



for i = 1:length(run)
    run(i).S.lon = [];
    run(i).S.lat = [];
    run(i).S.file_number = [];

    DIR = run(i).dir;
    D = dir(fullfile(DIR, '*high*.000*'));
    
    files = [];
    for year = years
        files = [files; dir(fullfile(run(i).dir, sprintf('*high*%04i*.nc',year)))];
    end

    for month = 1:length(files)
        fi = fullfile(files(month).folder, files(month).name);
        
        [~, ~, lon, lat] = cross_stream_sections_max(run(i).mesh_fi, fi, ssh_contour, lon_vec, lat_vec, 100);

        run(i).S.lon = cat(1, run(i).S.lon, lon);
        run(i).S.lat = cat(1, run(i).S.lat, lat);
        
        run(i).S.file_number = cat(1, run(i).S.file_number, month*ones(size(lon)));
    end
end

%% Plotting
i = 1;
figure(300+i)
clf
set(gcf,'color','w','position',[50 354 888 585])

xrange = [-77 -67];
yrange = [34 42];

m_proj('lambert','long',xrange,'lat',yrange);

m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linestyle','none','box','fancy','tickdir','out','fontsize',FS)

hold on
for k = 1:max(run(i).S.file_number)
    inds = run(i).S.file_number == k;

    m_line(run(i).S.lon(inds), run(i).S.lat(inds), 'color',run(i).color)

end


for j = 1:length(transect)
    m_line(transect(j).lon, transect(j).lat, 'color',rgb('red'),'linewidth',2,'marker','o','markerfacecolor','w')
end

title(sprintf('%s SSH %.1f m Years %i-%i', run(i).name, ssh_contour, years(1), years(end)))

set(gca,'fontsize',FS)

% m_legend(handle, run(:).name)
%%
saveas(gcf, fullfile(save_dir, sprintf('streamlines_%s_v1.png',run(i).name)))





