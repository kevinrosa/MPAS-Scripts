% Kevin Rosa
% July 29, 2019
%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = 'EC60to30_G_case';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).color = rgb('blue');
i = i+1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).color = rgb('black');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).color = rgb('red');


%% Settings
region = 'cal';

% download range:
xrange = [-131 -120];
yrange = [31 43];

% which .mat files to read:
target_string = 'timeall_lon-143.0to-116.0_lat_30.0to43.0_highFrequency';

% which fields to read from .mat files:
FIELDS3D = {'ssh','kineticEnergyAtSurface','temperatureAtSurface'};

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

%% convert KE to speed 
for i = 1:length(run)
    run(i).speed = sqrt(2 * run(i).kineticEnergyAtSurface);
    
end

%% pick a time for snapshots
target_date = datenum(0,07,1);
for i = 1:length(run)
    dv = datevec(run(i).time);
    day_of_year = datenum(0,dv(:,2),dv(:,3));
    
    close_times = find( abs(day_of_year-target_date) < 6 );
    
    % pick a date 3/4 through the run
    run(i).snap_t_ind = close_times(ceil(length(close_times)*0.5));
end

%% calculate grid cell widths
var_res_ind = 2; 
transition_contour = 12;

i = var_res_ind;
[~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;

%% sst gradient
for i = 1:length(run)
    run(i).sstgrad = NaN(size(run(i).temperatureAtSurface));
    for t = 1:length(run(i).time)
        run(i).sstgrad(t,:,:) = gradient_magnitude(run(i).LON, run(i).LAT, squeeze(run(i).temperatureAtSurface(t,:,:)));  % K/m
    end
end

%%
PLOTS = {'ssh_mean','ssh_snap','ssh_varlog10','speed_snap','speed_mean','sst_snap','sst_mean','sstgrad_snap','sstgrad_mean'};

% dict for determining which field is needed for each plot
FIELDS = containers.Map;
FIELDS('ssh') = 'ssh';
FIELDS('speed') = 'speed';
FIELDS('sst') = 'temperatureAtSurface';
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
P.(p).crange = [-1 1] * 0.25;
dc = 0.02;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 2*P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = flipud(cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));
P.(p).linecolor = 'k';
P.(p).cbar = 'SSH snapshot (m)';

p = 'ssh_mean';
P.(p).crange = [-1 1] * 0.25;
dc = 0.02;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 2*P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = flipud(cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));
P.(p).linecolor = 'k';
P.(p).cbar = 'SSH mean (m)';

p = 'ssh_varlog10';
P.(p).crange = [-3 -2];
dc = 0.1;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2)+1;
P.(p).cmap = cbrewer('seq','YlOrRd', length(P.(p).bins)-1, 'pchip');
P.(p).cmap(1,:) = [1,1,1];  % lowest color white
P.(p).linecolor = 'k';
P.(p).cbar = 'log10 of SSH variance (m)';

p = 'speed_snap';
P.(p).crange = [0 0.4];
dc = 0.001;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = flipud(cbrewer('div','Spectral', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));
P.(p).linecolor = 'none';
P.(p).cbar = 'Surface speed snapshot (m/s)';

p = 'speed_mean';
P.(p).crange = [0 0.25];
dc = 0.001;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = flipud(cbrewer('div','Spectral', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));
P.(p).linecolor = 'none';
P.(p).cbar = 'Surface speed mean (m/s)';

p = 'sst_snap';
P.(p).crange = [13 20];
dc = 0.1;
P.(p).bins = (P.(p).crange(1)-4):dc:P.(p).crange(2)+2;
P.(p).cmap = flipud(cbrewer('div','Spectral', length(P.(p).bins)-1, 'pchip'));
P.(p).linecolor = 'none';
P.(p).cbar = 'Sea Surface Temperature (^o C)';

p = 'sst_mean';
P.(p).crange = [13 20];
dc = 0.5;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = flipud(cbrewer('div','Spectral', length(P.(p).bins)-1, 'pchip'));
P.(p).linecolor = 'k';
P.(p).cbar = 'SST summer (JAS) mean (^o C)';

p = 'sstgrad_snap';
P.(p).crange = [0 3.2]*1e-5;
dc = 0.01*1e-5;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = parula(length(P.(p).bins)-1);
P.(p).linecolor = 'none';
P.(p).cbar = 'SST gradient (K/m)';

p = 'sstgrad_mean';
P.(p).crange = [0 2]*1e-5;
dc = 0.1*1e-5;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = parula(length(P.(p).bins)-1);
P.(p).linecolor = 'none';
P.(p).cbar = 'SST gradient summer (JAS) mean (K/m)';

%% plotting
w_mmap = 1;

run_inds = 1:3;

version_code = 'v2';

for plot_name = {'ssh_mean','ssh_snap','ssh_varlog10','sst_snap','sstgrad_snap','sstgrad_mean'}
    %%
    p = plot_name{1};

    for i = run_inds
    figure
    set(gcf,'color','w','position',[313 496 633 456])

    FIELD = P.(p).field;
    
    str = strsplit(p,'_');
    TITLE = [];
    
    if strcmp(P.(p).operation, 'snap')
        data = squeeze(run(i).(FIELD)(run(i).snap_t_ind,:,:));
        TITLE = sprintf('%s', datestr(run(i).time(run(i).snap_t_ind),'mm-dd-yyyy'));
    elseif strcmp(P.(p).operation, 'mean')
        if contains(p,'sst')
            % only average across summer months
            dv = datevec(run(i).time);
            t_inds = dv(:,2) >= 7 & dv(:,2) <= 9;
            data = squeeze(mean(run(i).(FIELD)(t_inds,:,:), 1));
        else
            data = squeeze(mean(run(i).(FIELD), 1));
        end
    elseif strcmp(P.(p).operation, 'var')
        data = squeeze(var(run(i).(FIELD), 1));
    elseif strcmp(P.(p).operation, 'varlog10')
        data = log10(squeeze(var(run(i).(FIELD), 1)));
    elseif strcmp(P.(p).operation, 'const')
        data = run(i).(FIELD);
    end
    
    if w_mmap
        m_proj('lambert','long',xrange,'lat',yrange);
        hold on

        m_contourf(run(i).LON, run(i).LAT, data, P.(p).bins,'linecolor',P.(p).linecolor)
    
        % add transition region contour
        if i == var_res_ind
            m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.2*[1,1,1],'linewidth',4)
            m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',1.0*[1,1,1],'linewidth',2)
        end

        m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
        m_grid('box','fancy','tickdir','out','fontsize',12)
        
    else
        contourf(run(i).LON, run(i).LAT, data, P.(p).bins,'linecolor',P.(p).linecolor)
        hold on    
        % add transition region contour
        if i == var_res_ind
            contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)
        end
    end

    caxis(P.(p).crange)
    colormap(P.(p).cmap)

    % SAVE
    save_name = sprintf('figures/mapviews/map_%s_%s_%s_%s', region, p, run(i).short_name, version_code);
    set(gca,'color','none');
%     export_fig(save_name, '-pdf','-painters','-transparent','-nocrop')
    export_fig(save_name, '-m3','-painters','-transparent','-nocrop')
    
    end
    

    %% save just colorbar
    figure 
    set(gcf,'color','w','position',[560 496 726 430])
    colormap(P.(p).cmap)
    caxis(P.(p).crange)

    cb = colorbar;
    ylabel(cb, P.(p).cbar)
    set(cb,'fontsize',14)
    set(gca,'visible','off')
    
    save_name = sprintf('figures/mapviews/map_%s_%s_colorbar_%s', region, p, version_code);
%     export_fig(save_name, '-pdf','-painters','-transparent')
    export_fig(save_name, '-m3','-painters','-transparent')
   
end







    
    
