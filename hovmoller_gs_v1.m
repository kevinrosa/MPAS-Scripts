% Kevin Rosa
% July 8, 2019

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
xrange = [-80 -30];
yrange = [34 36];

% which .mat files to read:
target_string = 'timeall_lon-97.0to-50.0_lat_18.0to45.0_highFrequency';

% which fields to read from .mat files:
FIELDS3D = {'ssh'};

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
    
    % lowpass filter ssh_trend
    cutoff_period_years = 5;
    cutoff_period = cutoff_period_years * 364.25 * 24;  % in hours

    
    [ssh_lp, time_lp] = lowpassfilter(run(i).time, run(i).ssh_trend, cutoff_period);
    run(i).ssh_trend_lp = interp1(time_lp, ssh_lp, run(i).time);
    
    run(i).ssh = run(i).ssh_raw - repmat(run(i).ssh_trend,[1,size(run(i).LON)]);
end

%%
figure
line(run(i).time, run(i).ssh_trend)
line(run(i).time, run(i).ssh_trend_lp)


%% calculate grid cell widths
var_res_ind = 2; 
transition_contour = 12;

i = var_res_ind;
[~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;



%%
target_lat = 34.8;

cmap = flipud(cbrewer('div','RdBu',60,'pchip'));

crange = [-1 1]*0.3;

yrange = [0 datenum(14,1,1)];
xrange = [-76 -55];

ref_speed = 0.04;  % m/s

for i = 1:length(run)
    figure(75+i)
    clf
    set(gcf,'color','w','position',[414 126 567 826])
    
    [~,eta] = min(abs(run(i).LAT(1,:) - target_lat));
    xi = run(i).LON(:,1)>xrange(1) & run(i).LON(:,1)<xrange(2);
    
    dx = lonlat_to_dxdy(run(i).LON(1,eta), run(i).LAT(1,eta), run(i).LON(2,eta), run(i).LAT(2,eta));  % km
    x = cumsum(repmat(dx, [length(run(i).LON(xi,eta)), 1]));  % km
    
    ssh = run(i).ssh_raw - repmat(run(i).ssh_trend,[1,size(run(i).LON)]);  % time-mean at each point (removes large-scale structure)
    
    time_mean = repmat(mean(ssh(:,xi,eta),1,'omitnan'), [length(run(i).time),1]);
    time = run(i).time-run(i).time(1);

%     time_mean = repmat(mean(ssh(:,:,eta),1,'omitnan'), [length(run(i).time),1]);
%     pcolor(squeeze(run(i).LON(:,eta)), time, ssh(:,:,eta) - time_mean); shading flat
    pcolor(x, time, ssh(:,xi,eta) - time_mean); shading flat
    
    
    colormap(jet)
    xlabel('Distance from coast (km)')
    ylabel('Time (years)')
    
    datetick('y','yy')
    
    set(gca,'fontsize',14,'ydir','reverse','xdir','normal')
    
    caxis(crange)
    ylim(yrange)
    
    colormap(cmap)
    
    % draw speed line
    time_sec = time*24*3600; 
    position = x(end) - (time_sec * ref_speed * 1e-3);
%     line(position, time, 'linewidth',2,'color',0.3*[1,1,1]);
    line([position(1) position(end)], [time(1) time(end)], 'linewidth',3,'color',0.5*[1,1,1]);
    
    % draw transition line
    if i == var_res_ind
        [c.lon,c.lat] = contour_line(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour);
        transition_lon = interp1(c.lat, c.lon, target_lat);
        
        transition_x = interp1(run(i).LON(xi,eta), x, transition_lon);
        line(transition_x*[1,1], yrange, 'color','k','linewidth',3,'linestyle','--')
    end
    
    
    save_name = sprintf('figures/westward_propogation/hovmoeller_atl_%s_v2.png',run(i).short_name);
    set(gca,'color','none');
    export_fig(save_name, '-m4','-painters','-transparent')
end

%% save a fig with colorbar
cb = colorbar;
ylabel(cb,'SSH anomaly (m)','fontsize',14)

save_name = 'figures/westward_propogation/hovmoeller_atl_colorbar_v2.png';
export_fig(save_name, '-m4','-painters','-transparent')