% RMS of surface velocities; gulf stream pathlines
% 
% Kevin Rosa
% July 29, 2019

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
% download range:
xrange = [-84 -54];
yrange = [30 44];

% which .mat files to read:
target_string = 'timeall_lon-97.0to-50.0_lat_18.0to45.0_highFrequency';

% which fields to read from .mat files:
FIELDS3D = {'ssh','kineticEnergyAtSurface'};


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
        run(i).ssh_trend(t) = mean(data(:),'omitnan');
    end
     
    run(i).ssh = run(i).ssh_raw - repmat(run(i).ssh_trend,[1,size(run(i).LON)]);
end

%% convert KE to speed and calculate rms
for i = 1:length(run)
    run(i).speed = sqrt(2 * run(i).kineticEnergyAtSurface);
    
    run(i).speed_rms = squeeze(rms(run(i).speed, 1));

end

%% calculate grid cell widths
var_res_ind = 2; 
transition_contour = 12;

i = var_res_ind;
[~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;


%% plot speed RMS
rms_cont= 0.4;

xrange = [-81 -56];
yrange = [32 43];

dc = 0.01;
crange = [0.3, 1.3];
conts1 = crange(1):dc:crange(2);
inds = conts1 < rms_cont;
cmap(inds,:) = repmat([1,1,1],length(find(inds)),1);
cmap(~inds,:) = flipud(cbrewer('div','Spectral',(length(find(~inds))), 'pchip'));
bins = 0:dc:2;

for i = 2:3
    figure
    set(gcf,'color','w','position',[560 565 695 361])

    m_proj('lambert','long',xrange,'lat',yrange);
    hold on
    m_contourf(run(i).LON, run(i).LAT, run(i).speed_rms, bins, 'linestyle','none')

    m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
    m_grid('box','fancy','tickdir','out','fontsize',12)
    
    if i == var_res_ind
        m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1], 'linecolor',0.6*[1,1,1],'linewidth',2)
    end

    colormap(cmap)
    caxis(crange)
    
    save_name = sprintf('figures/gulfstream/gs_rms_%s_v2.eps', run(i).short_name);
    set(gca,'color','none');
    export_fig(save_name,'-painters','-transparent','-nocrop')
end


%% save just colorbar
figure 
set(gcf,'color','w')
colormap(cmap)
caxis(crange)

cb = colorbar;
ylabel(cb, 'Surface speed RMS (m/s)')
set(cb,'fontsize',14)
cbarrow
set(gca,'visible','off')

save_name = 'figures/gulfstream/gs_rms_colorbar_v2.eps';
export_fig(save_name,'-painters','-transparent')


%%
%
%
%
%
%% Gulf stream pathlines
% run(i).S will hold streamline info

no_years = 3;
ssh_contour = -0.2;
run_inds = 2:3;  % run indices to run calculation on

for i = run_inds
    
    time_range = [run(i).time(end)-datenum(no_years,1,1), run(i).time(end)];
    t_inds = find(run(i).time >= time_range(1) & run(i).time <= time_range(end));
    
    t = 1;  % gulf stream contour counter
    for t_ind = 1:6:length(t_inds)  % only plot every 30 days
        run_t = t_inds(t_ind);
        
        [lon, lat] = streamline_coords(run(i).LON, run(i).LAT, squeeze(run(i).ssh(t,:,:)), ssh_contour);
        run(i).gs(t).lon = lon;
        run(i).gs(t).lat = lat;
        
        t = t+1;
    end
    
end


%% Plotting streamlines
figure
set(gcf,'color','w','position',[560 565 695 361])

m_proj('lambert','long',xrange,'lat',yrange);
hold on

for i = [3,2]
    for t = 1:length(run(i).gs)
         m_line(run(i).gs(t).lon, run(i).gs(t).lat, 'color',run(i).color)
    end
end

i = var_res_ind;
m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1], 'linecolor',0.6*[1,1,1],'linewidth',4)

m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('box','fancy','tickdir','out','fontsize',12)


save_name = 'figures/gulfstream/gs_paths_v2.eps';
set(gca,'color','none');
export_fig(save_name, '-painters','-transparent','-nocrop')



