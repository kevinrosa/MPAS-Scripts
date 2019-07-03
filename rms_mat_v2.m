% RMS of surface vvelocities
% 
% Kevin Rosa
% May 29, 2019

addpath(genpath('.'))

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 15:20;%6:22;
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).years = 25:30;%20:36;
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = '20180305.GM600.T62_oECv3.eos';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).years = 25:30;%20:36;
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
xrange = [-97 -74];
yrange = [18 35];
dx = 0.1;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

FIELDS = {'speed','ssh','temperatureAtSurface','temperatureAt250m','relativeVorticityAt250m'};%,'dThreshMLD','tThreshMLD'};

for i = 1:length(run)
    files = {};
    for year = run(i).years
        dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.highFrequencyOutput.%04i-*',year)));
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
        
        run(i).time(tt) = mpas_time(data_fi, t_ind);            

        
        % Calculate LON LAT and land-sea mask (only do this once)
        if tt == 1
            [LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
            run(i).LON = LON;
            run(i).LAT = LAT;
            
            run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);
        end
        
        % read data fields
        for F = FIELDS
            if strcmp('speed',F{1})
                [~,~,kineticenergy] = mpas_to_lonlat_meshgrid('kineticEnergyAtSurface', run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
                data = sqrt(2 * kineticenergy);
            else
                [~,~,data] = mpas_to_lonlat_meshgrid(F{1}, run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
            end
            run(i).(F{1})(tt,:,:) = data .* run(i).mask;
        end
        
%         ref_ssh = interp1(run(i).S.time, run(i).S.ssh_mean, run(i).time(tt));
%         run(i).ssh(tt,:,:)  = run(i).ssh(tt,:,:) - ref_ssh;
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
        tt = tt+1;
    end
    
    run(i).speed_rms = squeeze(rms(run(i).speed, 1));
end

%% De-trend SSH
for i = 1:length(run)
    ss = reshape(run(i).ssh, length(run(i).time), []);
    run(i).ssh_mean = mean(ss,2,'omitnan');
    run(i).SSH = run(i).ssh - repmat(run(i).ssh_mean, [1, size(run(i).LON)]);
end

%%
save_dir = 'figures/rms_mat_v2';

%% 2d field (no time dimension
rms_cont= 0.4;
for fields = {'speed_rms'}
figure
clf; set(gcf,'color','w','position',[41 236 1576 697])
field = fields{1};

run_inds = [3,1,2];

if strcmp('temperatureAtSurface', field)
    crange = [10 20];
elseif strcmp('relativeVorticityAt250m',field)
    crange = 1e-5 * [-1 1];
elseif strcmp('ssh',field)
    dc = 0.015;
    bins = -0.3:dc:0;
    crange = [bins(1) bins(end)];
elseif strcmp('speed_rms',field)
    crange = [0 1.5];
end

m = 1;
n = length(run_inds);
p = 1;  % subplot counter
for i = run_inds   
    
    subplot(m,n,p)
    
    m_proj('lambert','long',xrange,'lat',yrange);
    m_pcolor(run(i).LON, run(i).LAT, run(i).speed_rms)
    hold on
    m_contour(run(i).LON, run(i).LAT, run(i).speed_rms, rms_cont*[1,1], 'color','w','linewidth',2)

    m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
    m_grid('linestyle','none','box','fancy','tickdir','out','fontsize',12)

    caxis(crange)
    colormap(jet)
    cb = colorbar;
    ylabel(cb,'Surface Speed RMS (m/s)')
        
    title(sprintf('%s %s', field, run(i).short_name),'interpreter','none')
    p = p+1;
end

save_name = fullfile(save_dir, sprintf('3subplots_%s_cont%.2f_001.png',field, rms_cont));
saveas(gcf, save_name)
end

%% mean
rms_cont= 0.4;
for fields = {'speed'}
figure
clf; set(gcf,'color','w','position',[41 236 1576 697])
field = fields{1};

run_inds = [3,1,2];

if strcmp('temperatureAtSurface', field)
    crange = [10 20];
elseif strcmp('relativeVorticityAt250m',field)
    crange = 1e-5 * [-1 1];
elseif strcmp('ssh',lower(field))
    dc = 0.025;
    crange = 0.6 * [-1 1];
    bins = crange(1):dc:crange(2);
elseif strcmp('speed_rms',field)
    crange = [0 1.5];
elseif strcmp('speed',field)
    crange = [0 1.5];    
end

m = 1;
n = length(run_inds);
p = 1;  % subplot counter
for i = run_inds   
    
    subplot(m,n,p)
    
    m_proj('lambert','long',xrange,'lat',yrange);
    m_pcolor(run(i).LON, run(i).LAT, squeeze(mean(run(i).(field),1)))
    hold on
%     m_contour(run(i).LON, run(i).LAT, squeeze(mean(run(i).(field),1)), bins, 'linecolor','k')
%     hold on
    
    if contains(field,'speed')
        m_contour(run(i).LON, run(i).LAT, squeeze(mean(run(i).(field),1)), rms_cont*[1,1], 'color','w','linewidth',2)
        colormap(jet)
    end

    m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
    m_grid('linestyle','none','box','fancy','tickdir','out','fontsize',12)

    caxis(crange)
    cb = colorbar;
    ylabel(cb,'Mean Surface Speed (m/s)')
        
    title(sprintf('Mean %s %s', field, run(i).short_name),'interpreter','none')
    p = p+1;
end

if contains(field,'speed')
    save_name = fullfile(save_dir, sprintf('3subplots_mean_%s_cont%.2f_001.png',field, rms_cont));
else
    save_name = fullfile(save_dir, sprintf('3subplots_mean_%s_cont%.3f_002.png',field, dc));
end
saveas(gcf, save_name)
end

%%
%% Download Grid resolution
for i = 1:length(run)
    [~,~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, lon_vec, lat_vec, 0);
    
    run(i).dx = 2 * sqrt(run(i).mask .* areaCell ./ pi) * 1e-3;  % (km)
    
end

%% Plot resolution
figure
clf; set(gcf,'color','w')%,'position',[121 357 1363 594])

run_inds = 1;%[3,1,2];

dc = 5;
bins = 5:dc:70;

m = 1;
n = length(run_inds);
p = 1;  % subplot counter
for i = run_inds   
    
    subplot(m,n,p)
    contourf(run(i).LON, run(i).LAT, run(i).dx, bins);
    colorbar
    
%     caxis(crange)
    
    daspect([1, cosd(mean(run(i).LAT(:))), 1])
    
    title(sprintf('Grid cell width (km) %s', run(i).short_name))
    p = p+1;
end

save_name = sprintf('figures/cali_001/grid_resolution_%s_cont%ikm_001.png',run(i).short_name,dc);
saveas(gcf, save_name)











%%
xrange = [-85 -55];
yrange = [lat_vec(1), lat_vec(end)];

FS = 14;
crange = [0 1.5];

rms_cont = 0.4;
white_cont = 0.4;

%% convoluted way to make colormap 
conts = crange(1):0.05:crange(end);
cmap = jet(length(conts)-1);
% vals = (conts(1:end-1)+conts(2:end))/2;
% inds = vals<=white_cont;
% cmap(inds,:) = repmat([1,1,1],length(find(inds)),1);
% cmap(~inds,:) = jet(length(find(~inds)));

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

colormap(cmap)

saveas(gcf, fullfile(save_dir, sprintf('speed_rms_pcolor_%s_showtransects_v0.png', run(i).short_name)))
end

%%
figure(210)
clf
set(gcf,'color','w','position',[50 354 888 585])
m_proj('lambert','long',xrange,'lat',yrange);

m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linestyle','none','box','fancy','tickdir','out','fontsize',FS)

hold on
for i = 1:2%length(run)
    
    [~,handle(i)] = m_contour(run(i).lon, run(i).lat, run(i).speed_rms, rms_cont*[1,1], 'color',run(i).color,'linewidth',2);
    
end

for j = 1:2%length(transect)
    m_line(transect(j).lon, transect(j).lat, 'color',rgb('red'),'linewidth',2,'marker','o','markerfacecolor','w')
end

title(sprintf('%.1f m/s RMS Contours: %s vs. %s',rms_cont,run(1).short_name,run(2).short_name))

set(gca,'fontsize',FS)

% m_legend(handle, run(:).name)

saveas(gcf, fullfile(save_dir, sprintf('speed_rms_contours_showtransects_v1.png')))

%%
%
%
%
%
%% Gulf stream pathlines
% run(i).S will hold streamline info

years = 7:14;
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
for i = 1:length(run)
    figure(300+i)
    clf
    set(gcf,'color','w','position',[50 354 888 585])

    xrange = [-77 -64];%-67];
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

    pause(0.1)
    title(sprintf('%s SSH %.1f m Years %i-%i', run(i).name, ssh_contour, years(1), years(end)))

    set(gca,'fontsize',FS)
    
    saveas(gcf, fullfile(save_dir, sprintf('streamlines_%s_v0.png',run(i).short_name)))
end





