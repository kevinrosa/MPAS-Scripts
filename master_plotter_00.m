% Kevin Rosa
% June 12, 2019
%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = '20180305.GM600.T62_oECv3.eos';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).color = rgb('blue');

% matObj = matfile('Pacific_years20to36_time1_lon-143.0to-116.0_lat_30.0to43.0_highFrequency_20180208.GMPAS-IAF.T62_oRRS18v3.anvil.mat');
% varlist = who(matObj)

%% Settings
region = 'cc0';

% download range:
xrange = [-131 -120];
yrange = [31 43];

% which .mat files to read:
target_string = 'time1_lon-143.0to-116.0_lat_30.0to43.0_highFrequency';

% which fields to read from .mat files:
% FIELDS3D = {'ssh','temperatureAtSurface','barotropicSpeed','relativeVorticityAt250m'};
FIELDS3D = {'ssh','temperatureAtSurface'};

%% Load
for i = 1:length(run)
    D = dir(sprintf('*%s*%s.mat', target_string, run(i).code));
    run(i).fi = D.name;


    % determine spatial indices
    tmp = load(run(i).fi, 'LON','LAT');
    xi  = tmp.LON(:,1)>=xrange(1) & tmp.LON(:,1) <=xrange(2);
    eta = tmp.LAT(1,:)>=yrange(1) & tmp.LAT(1,:) <=yrange(2);

    run(i).LON = tmp.LON(xi,eta);
    run(i).LAT = tmp.LAT(xi,eta);

    % 1D fields
    for F = {'time'}
        tmp = load(run(i).fi, F{1});
        run(i).(F{1}) = tmp.(F{1});
    end

    % additional 2D fields
    for F = {'mask'}
        tmp = load(run(i).fi, F{1});
        run(i).(F{1}) = tmp.(F{1})(xi,eta);
    end

    % 3D fields
    for F = FIELDS3D
        tmp = load(run(i).fi, F{1});
        run(i).(F{1}) = tmp.(F{1})(:,xi,eta);
    end
    
    % pre de-trended ssh
    run(i).ssh_raw = run(i).ssh;
end

%% de-trend SSH
for i = 1:length(run)
    run(i).ssh_trend = NaN(length(run(i).time),1);
    for t = 1:length(run(i).time)
        data = run(i).ssh_raw(t,:,:);
        run(i).ssh_trend(t) = nanmean(data(:));
    end
    
    run(i).ssh = run(i).ssh_raw - repmat(run(i).ssh_trend,[1,size(run(i).LON)]);
end

%% pick a time for snapshots
target_date = datenum(0,07,1);
for i = 1:length(run)
    dv = datevec(run(i).time);
    day_of_year = datenum(0,dv(:,2),dv(:,3));
    
    close_times = find( abs(day_of_year-target_date) < 6 );
    
    % pick a date in the middle of the run
    run(i).snap_t_ind = close_times(ceil(length(close_times)/2));
end

%% calculate grid cell widths
for i = 1:length(run)
    [~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
    run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;
end

%% calculate meridional geostrophic velocities
for i = 1:length(run)
    run(i).vgeo = NaN(size(run(i).ssh));
    
    for t = 1:length(run(i).time)
        run(i).vgeo(t,:,:) = geostrophic_vel_meridional(run(i).LON, run(i).LAT, squeeze(run(i).ssh(t,:,:)));
    end
end

%% sst gradients
for i = 1:length(run)
    run(i).sstgrad = NaN(size(run(i).temperatureAtSurface));
    
    for t = 1:length(run(i).time)
        run(i).sstgrad(t,:,:) = gradient_magnitude(run(i).LON, run(i).LAT, squeeze(run(i).temperatureAtSurface(t,:,:)));
    end
end

%%
PLOTS = {'sst_grad'};%{'res_const','ssh_mean','ssh_var','ssh_snap','sst_mean','sst_snap','sst_var','btspeed_mean','btspeed_snap','btspeed_var','vort_snap','vort_mean','vgeo_snap','vgeo_mean','vgeo_var'};
PLOTS = {'sstgrad_mean'};

% dict for determining which field is needed for each plot
FIELDS = containers.Map;
FIELDS('ssh') = 'ssh';
FIELDS('sst') = 'temperatureAtSurface';
FIELDS('btspeed') = 'barotropicSpeed';
FIELDS('vort') = 'relativeVorticityAt250m';
FIELDS('res') = 'widthCell';
FIELDS('vgeo') = 'vgeo';
FIELDS('sstgrad') = 'sstgrad';

for plot_name = PLOTS
    p = plot_name{1};
    str = strsplit(p ,'_');
    P.(p).field = FIELDS(str{1});
end

% determine which operation (mean, variance, snapshot, etc.)
for plot_name = PLOTS
    p = plot_name{1};
    str = strsplit(p ,'_');
    P.(p).operation = str{end};
end


%%
p = 'ssh_snap';
P.(p).crange = [-1 1] * 0.3;
dc = 0.01;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = flipud(cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));

p = 'ssh_mean';
P.(p).crange = [-1 1] * 0.25;
dc = 0.01;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = flipud(cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));

p = 'ssh_var';
P.(p).crange = [0 1] * 0.01;
dc = 0.0005;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cbrewer('seq','YlOrRd', length(P.(p).bins)-1, 'pchip');
P.(p).cmap(1,:) = [1,1,1];

p = 'sst_snap';
P.(p).crange = [13 20];
dc = 0.1;
P.(p).bins = (P.(p).crange(1)-4):dc:P.(p).crange(2);
P.(p).cmap = cmocean('thermal',length(P.(p).bins)-1);
P.(p).cmap = flipud(cbrewer('div','Spectral', length(P.(p).bins)-1, 'pchip'));

p = 'sst_mean';
P.(p).crange = [10 20];
dc = 0.5;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cmocean('thermal',length(P.(p).bins)-1);

p = 'sst_var';
P.(p).crange = [0 10];
dc = 0.25;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cbrewer('seq','YlOrRd', length(P.(p).bins)-1, 'pchip');
P.(p).cmap(1,:) = [1,1,1];

p = 'btspeed_mean';
P.(p).crange = [0 0.1];
dc = 0.001;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = flipud(cbrewer('div','Spectral', length(P.(p).bins)-1, 'pchip'));

p = 'btspeed_snap';
P.(p).crange = [0 0.1];
dc = 0.001;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = flipud(cbrewer('div','Spectral', length(P.(p).bins)-1, 'pchip'));

p = 'btspeed_var';
P.(p).crange = [0 0.0004];
dc = 0.00001;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cbrewer('seq','YlOrRd', length(P.(p).bins)-1, 'pchip');
P.(p).cmap(1,:) = [1,1,1];

p = 'vort_snap';
P.(p).crange = [-1 1]*5e-6;
dc = 1e-7;
P.(p).bins = 4*P.(p).crange(1):dc:P.(p).crange(2)*4;
P.(p).cmap = flipud(cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));

p = 'vort_mean';
P.(p).crange = [-1 1]*1e-6;
dc = 1e-7;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 2*P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = flipud(cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));

p = 'vort_var';
P.(p).crange = [0 1]*1e-4;
dc = P.(p).crange(2) / 20;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cbrewer('seq','YlOrRd', length(P.(p).bins)-1, 'pchip');
P.(p).cmap(1,:) = [1,1,1];

p = 'res_const';
P.(p).crange = [6 60];
dc = 1;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cbrewer('div','Spectral', length(P.(p).bins)-1, 'pchip');

p = 'vgeo_snap';
P.(p).crange = [-1 1]*0.25;
dc = 0.01;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 5*P.(p).crange(1):dc:P.(p).crange(2)*5 ;
P.(p).cmap = cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip');

p = 'vgeo_mean';
P.(p).crange = [-1 1]*0.1;
dc = 0.02;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 2*P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip');

p = 'vgeo_var';
P.(p).crange = [0 0.015];
dc = 0.0005;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cbrewer('seq','YlOrRd', length(P.(p).bins)-1, 'pchip');
P.(p).cmap(1,:) = [1,1,1];

p = 'sst_grad';
P.(p).crange = [0 3.5]*1e-5;
dc = 0.1*1e-5;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = parula(length(P.(p).bins)-1);

p = 'sstgrad_mean';
P.(p).crange = [0 1]*1e-5;
dc = 0.1*1e-5;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = parula(length(P.(p).bins)-1);

%% w/ m_map
run_inds = [3,1,2];
m = 1;
n = length(run_inds);

version_code = 'v0';

var_res_ind = 1;
transition_contour = 15;

for plot_name = {'sstgrad_mean'}%{'sst_grad'}%{'sst_snap'}%{'vort_snap','vort_mean','ssh_snap','ssh_mean'}%{'vgeo_snap','vgeo_mean','vgeo_var'}% PLOTS
    p = plot_name{1};
    
    figure
    set(gcf,'name',p,'position',[13 449 1629 503],'color','w')
    pcounter = 1;

    for i = run_inds
    subplot(m,n,pcounter)

    FIELD = P.(p).field;
    
    str = strsplit(p,'_');
    TITLE = sprintf('%s %s %s', run(i).short_name, str{1}, P.(p).operation);
    
    if strcmp(P.(p).operation, 'snap')
        data = squeeze(run(i).(FIELD)(run(i).snap_t_ind,:,:));
        TITLE = sprintf('%s %s', TITLE, datestr(run(i).time(run(i).snap_t_ind),'mm-dd-yyyy'));
    elseif strcmp(P.(p).operation, 'mean')
        data = squeeze(mean(run(i).(FIELD), 1));
    elseif strcmp(P.(p).operation, 'var')
        data = squeeze(var(run(i).(FIELD), 1));
    elseif strcmp(P.(p).operation, 'const')
        data = run(i).(FIELD);
    elseif strcmp(P.(p).operation, 'grad')
        data = gradient_magnitude(run(i).LON, run(i).LAT, squeeze(run(i).(FIELD)(run(i).snap_t_ind,:,:)));
        TITLE = sprintf('%s %s', TITLE, datestr(run(i).time(run(i).snap_t_ind),'mm-dd-yyyy'));
    end
    
    m_proj('lambert','long',xrange,'lat',yrange);
    hold on
    
    if strcmp(p,'btspeed_mean') | strcmp(p,'btspeed_snap') | strcmp(p,'vort_snap') | strcmp(p,'vgeo_snap') | strcmp(p,'sst_snap') | strcmp(p,'sst_grad') | strcmp(p,'sstgrad_mean')
        m_contourf(run(i).LON, run(i).LAT, data, P.(p).bins,'linecolor','none')
    else
        m_contourf(run(i).LON, run(i).LAT, data, P.(p).bins,'linecolor','k')
    end
    
%     cont_level = 16.2;
%     if strcmp(p,'sst_snap')
%         m_contour(run(i).LON, run(i).LAT, data, cont_level*[1,1],'linecolor','k')
%     end
    
    % add transition region contour
    if i == var_res_ind
        m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)
    end
    
    m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
    m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')

    caxis(P.(p).crange)
    colormap(P.(p).cmap)
    colorbar
    
    title(TITLE)

    pcounter = pcounter+1;
    end
    
    save_name = sprintf('figures/master_plotter_00/%s_%s_%s.png', region, p, version_code);
%     save_name = sprintf('figures/master_plotter_00/%s_%s_cont%.1fspectral_%s.png', region, p, cont_level, version_code);    
    saveas(gcf, save_name)
end


 