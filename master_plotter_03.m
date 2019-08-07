% Kevin Rosa
% July 17, 2019

%%
addpath(genpath('.'))

%%
type = 'mpaso.hist.am.highFrequencyOutput';

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).color = rgb('red');
run(i).years = 20;
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).color = rgb('black');
run(i).years = 20;
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = 'EC60to30_G_case';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).color = rgb('blue');
run(i).years = 50;

% matObj = matfile('Pacific_years20to36_time1_lon-143.0to-116.0_lat_30.0to43.0_highFrequency_20180208.GMPAS-IAF.T62_oRRS18v3.anvil.mat');
% varlist = who(matObj)

%% Settings
region = 'na0';

% download range:
xrange = [-80 20];
yrange = [45 80];

%%
FIELDS = {'dThreshMLD','tThreshMLD'};
t_ind = 'all';
months = 1:3;

%%
dx = 0.1;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);


for i = 1:length(run)
    files = {};
    for year = run(i).years
        for month = months
            dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.highFrequencyOutput.%04i-%02i-*',year,month)));
            files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
        end
    end

%     t_inds_per_year = 350/5;  % underestimate a few since it's better to be under than over
%     no_years = length(files)/12;
    t_length = 10; %floor(t_inds_per_year * no_years);  % approx number of time indices -- used to initialize arrays 

    nans1d = NaN(t_length, 1);
    run(i).time = nans1d;
    
    nans3d = NaN([t_length, length(lon_vec), length(lat_vec)]);
    for F = FIELDS
        run(i).(F{1}) = nans3d;
    end

    t1 = 1;  % time index counter
    for m = 1:length(files)
        
        data_fi = files{m};
        
        new_time = mpas_time(data_fi, t_ind);
        t2 = t1 + length(new_time) - 1;
        
        run(i).time(t1:t2) = new_time;            
        
        % Calculate LON LAT and land-sea mask (only do this once)
        if t1 == 1
            [LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
            run(i).LON = LON;
            run(i).LAT = LAT;
            
            run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);
        end
        
        % read data fields
        for F = FIELDS
            [~,~,data] = mpas_to_lonlat_meshgrid(F{1}, run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);        
            run(i).(F{1})(t1:t2,:,:) = data;

        end
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*m/length(files), '%')
        t1 = t2+1;
    end
    
    % apply mask
    for F = FIELDS
        mask = repmat(permute(run(i).mask,[3,1,2]), [length(run(i).time),1,1]);
        run(i).(F{1}) = run(i).(F{1}) .* mask;
    end
    
%     % save data to mat
%     save_name = sprintf('%s_years%ito%i_time%s_lon%.1fto%.1f_lat_%.1fto%.1f_%s_%s.mat',...
%         region(j).name, run(i).years(1), run(i).years(end), num2str(t_ind),...
%         region(j).lon_range(1), region(j).lon_range(2), region(j).lat_range(1), region(j).lat_range(2),...
%         type, run(i).code);
%     out = run(i);
%     
%     save(save_name, '-struct', 'out');  % -struct saves structure fields as separate variables in the mat file 
%     
%     % reduce memory usage a little:
%     clear out
%     for F = FIELDS
%         run(i).(F{1}) = [];
%     end
end


%%
% 
% % which .mat files to read:
% target_string = 'time1_lon-97.0to-50.0_lat_18.0to45.0_highFrequency';
% 
% % which fields to read from .mat files:
% FIELDS3D = {'ssh','temperatureAtSurface','barotropicSpeed','relativeVorticityAt250m'};
% 
% %% Load
% for i = 1:length(run)
%     D = dir(sprintf('*%s*%s.mat', target_string, run(i).code));
%     run(i).fi = D.name;
% 
% 
%     % determine spatial indices
%     tmp = load(run(i).fi, 'LON','LAT');
%     xi  = tmp.LON(:,1)>=xrange(1) & tmp.LON(:,1) <=xrange(2);
%     eta = tmp.LAT(1,:)>=yrange(1) & tmp.LAT(1,:) <=yrange(2);
% 
%     run(i).LON = tmp.LON(xi,eta);
%     run(i).LAT = tmp.LAT(xi,eta);
% 
%     % 1D fields
%     for F = {'time'}
%         tmp = load(run(i).fi, F{1});
%         run(i).(F{1}) = tmp.(F{1});
%     end
% 
%     % additional 2D fields
%     for F = {'mask'}
%         tmp = load(run(i).fi, F{1});
%         run(i).(F{1}) = tmp.(F{1})(xi,eta);
%     end
% 
%     % 3D fields
%     for F = FIELDS3D
%         tmp = load(run(i).fi, F{1});
%         run(i).(F{1}) = tmp.(F{1})(:,xi,eta);
%     end
%     
%     % pre de-trended ssh
%     run(i).ssh_raw = run(i).ssh;
% % end
% 
% %% de-trend SSH
% for i = 1:length(run)
%     run(i).ssh_trend = NaN(length(run(i).time),1);
%     for t = 1:length(run(i).time)
%         data = run(i).ssh_raw(t,:,:);
%         run(i).ssh_trend(t) = nanmean(data(:));
%     end
%     
%     run(i).ssh = run(i).ssh_raw - repmat(run(i).ssh_trend,[1,size(run(i).LON)]);
% end

% %% pick a time for snapshots
% target_date = datenum(0,07,1);
% for i = 1:length(run)
%     dv = datevec(run(i).time);
%     day_of_year = datenum(0,dv(:,2),dv(:,3));
%     
%     close_times = find( abs(day_of_year-target_date) < 6 );
%     
%     % pick a date in the middle of the run
%     run(i).snap_t_ind = close_times(ceil(length(close_times)/2));
% end

%% calculate grid cell widths
for i = 1:length(run)
    [~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
    run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;
end

var_res_ind = 1; 
transition_contour = 15;

%% plot
FIELD = 'dThreshMLD';

figure
set(gcf,'name',FIELD,'position',[13 449 1629 503],'color','w')

run_inds = [3,1,2];
m = 1;
n = 3;
pcounter = 1;
for i = run_inds
subplot(m,n,pcounter)

data = squeeze(mean(run(i).(FIELD), 1));


m_proj('lambert','long',xrange,'lat',yrange);
hold on

% m_contourf(run(i).LON, run(i).LAT, data, bins,'linecolor','none')
m_pcolor(run(i).LON, run(i).LAT, data); shading flat
% pcolor(run(i).LON, run(i).LAT, data); shading flat

hold on
% add transition region contour
if i == var_res_ind
    m_contour(run(i).LON, run(i).LAT, run(i).widthCell .* run(i).mask, transition_contour*[1,1],'linecolor',0.8*[1,1,1],'linewidth',2)
%     contour(run(i).LON, run(i).LAT, run(i).widthCell .* run(i).mask, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)
end

m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')


% caxis(P.(p).crange)
% colormap(P.(p).cmap)
caxis([0 1000])
colorbar

% title(TITLE)

pcounter = pcounter+1;
end

save_name = sprintf('figures/mld/northatlantic_mld_v1.png');%, region, p, version_code);
export_fig(save_name, '-m2')

%%
%
%%
PLOTS = {'res_const','ssh_mean','ssh_var','ssh_snap','sst_mean','sst_snap','sst_var','btspeed_mean','btspeed_snap','btspeed_var','vort_snap','vort_mean','vgeo_snap','vgeo_mean','vgeo_var'};


% dict for determining which field is needed for each plot
FIELDS = containers.Map;
FIELDS('ssh') = 'ssh';
FIELDS('sst') = 'temperatureAtSurface';
FIELDS('btspeed') = 'barotropicSpeed';
FIELDS('vort') = 'relativeVorticityAt250m';
FIELDS('res') = 'widthCell';
FIELDS('vgeo') = 'vgeo';

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
P.(p).crange = [-1 1] * 0.5;
dc = 0.02;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 5*P.(p).crange(1):dc:P.(p).crange(2)*5;
P.(p).cmap = flipud(cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));

p = 'ssh_mean';
P.(p).crange = [-1 1] * 0.5;
dc = 0.02;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 5*P.(p).crange(1):dc:P.(p).crange(2)*5;
P.(p).cmap = flipud(cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));

p = 'ssh_var';
P.(p).crange = [0 1] * 0.05;
dc = 0.001;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cbrewer('seq','YlOrRd', length(P.(p).bins)-1, 'pchip');
P.(p).cmap(1,:) = [1,1,1];

p = 'sst_snap';
P.(p).crange = [20 30];
dc = 0.15;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cmocean('thermal',length(P.(p).bins)-1);

p = 'sst_mean';
P.(p).crange = [20 30];
dc = 0.5;
P.(p).bins = 10:dc:P.(p).crange(2);
P.(p).cmap = cmocean('thermal',length(P.(p).bins)-1);

p = 'sst_var';
P.(p).crange = [0 15];
dc = 0.5;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cbrewer('seq','YlOrRd', length(P.(p).bins)-1, 'pchip');
P.(p).cmap(1,:) = [1,1,1];

p = 'btspeed_mean';
P.(p).crange = [0 1];
dc = 0.001;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = flipud(cbrewer('div','Spectral', length(P.(p).bins)-1, 'pchip'));

p = 'btspeed_snap';
P.(p).crange = [0 1];
dc = 0.001;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = flipud(cbrewer('div','Spectral', length(P.(p).bins)-1, 'pchip'));

p = 'btspeed_var';
P.(p).crange = [0 0.01];
dc = 0.0005;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cbrewer('seq','YlOrRd', length(P.(p).bins)-1, 'pchip');
P.(p).cmap(1,:) = [1,1,1];

p = 'vort_snap';
P.(p).crange = [-1 1]*1e-5;
dc = 1e-7;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 5*P.(p).crange(1):dc:P.(p).crange(2)*5;
P.(p).cmap = flipud(cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip'));

p = 'vort_mean';
P.(p).crange = [-1 1]*1e-6;
dc = 1e-7;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 5*P.(p).crange(1):dc:P.(p).crange(2)*5;
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
P.(p).crange = [-1 1]*0.2;
dc = 0.01;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 2*P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip');

p = 'vgeo_mean';
P.(p).crange = [-1 1]*0.1;
dc = 0.01;
P.(p).crange = [P.(p).crange(1)-dc/2, P.(p).crange(2)+dc/2];
P.(p).bins = 2*P.(p).crange(1):dc:P.(p).crange(2)*2;
P.(p).cmap = cbrewer('div','RdBu', length(P.(p).crange(1):dc:P.(p).crange(2))-1, 'pchip');

p = 'vgeo_var';
P.(p).crange = [0 0.01];
dc = 0.0005;
P.(p).bins = P.(p).crange(1):dc:P.(p).crange(2);
P.(p).cmap = cbrewer('seq','YlOrRd', length(P.(p).bins)-1, 'pchip');
P.(p).cmap(1,:) = [1,1,1];

%% plotting
w_mmap = 1;

run_inds = [3,1,2];
m = 1;
n = length(run_inds);

version_code = 'v0';

var_res_ind = 1;
transition_contour = 15;

for plot_name = {'vort_mean'}%,'vort_snap','vgeo_snap','vgeo_mean'}%,'ssh_snap','ssh_mean'}%PLOTS'vgeo_snap','vgeo_mean','vgeo_var',
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
    end
    
    if w_mmap
        m_proj('lambert','long',xrange,'lat',yrange);
        hold on

        if strcmp(p,'btspeed_mean') | strcmp(p,'btspeed_snap') | strcmp(p,'vort_snap') | strcmp(p,'vort_mean') | strcmp(p,'vgeo_snap')
            m_contourf(run(i).LON, run(i).LAT, data, P.(p).bins,'linecolor','none')
        else
            m_contourf(run(i).LON, run(i).LAT, data, P.(p).bins,'linecolor','k')
        end
    
        % add transition region contour
        if i == var_res_ind
            m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)
        end

        m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
        m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')
        
    else
        if strcmp(p,'btspeed_mean') | strcmp(p,'btspeed_snap') | strcmp(p,'vort_snap')
            contourf(run(i).LON, run(i).LAT, data, P.(p).bins,'linecolor','none')
        else
            contourf(run(i).LON, run(i).LAT, data, P.(p).bins,'linecolor','k')
        end
        
        hold on
    
        % add transition region contour
        if i == var_res_ind
            contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)
        end
    end

    caxis(P.(p).crange)
    colormap(P.(p).cmap)
    colorbar
    
    title(TITLE)

    pcounter = pcounter+1;
    end
    
    save_name = sprintf('figures/master_plotter_02/%s_%s_%s.png', region, p, version_code);
    saveas(gcf, save_name)
end


%% hovmoller
i = 2;
target_lat = 34.8;
[~,eta] = min(abs(run(i).LAT(1,:) - target_lat));

figure
field = 'temperatureAtSurface';
% data = run(i).(field)(:,:,eta) - repmat(mean(run(i).(field)(:,:,eta),1),[length(run(i).time),1,1]);

data = run(i).(field)(:,:,eta) - repmat(mean(run(i).(field)(:,:,eta),2,'omitnan'),[1,length(run(i).LON(:,1)),1]);


pcolor(run(i).LON(:,eta), run(i).time, data); shading flat

set(gca,'ydir','reverse')

colormap(parula)
%%
caxis([-0.5 0.5])
%%
caxis(4*[-1 1])




    
    
