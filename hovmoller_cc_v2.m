% Kevin Rosa
% July 31, 2019

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
xrange = [-143 -122];
yrange = [39 41];

% which .mat files to read:
target_string = 'timeall_lon-143.0to-116.0_lat_30.0to43.0_highFrequency';

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

    
%     [ssh_lp, time_lp] = lowpassfilter(run(i).time, run(i).ssh_trend, cutoff_period);
%     run(i).ssh_trend_lp = interp1(time_lp, ssh_lp, run(i).time);
%     
%     run(i).ssh = run(i).ssh_raw - repmat(run(i).ssh_trend,[1,size(run(i).LON)]);
end


%% calculate grid cell widths
var_res_ind = 2; 
transition_contour = 12;

i = var_res_ind;
[~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;



%%
target_lat = 40;

cmap = flipud(cbrewer('div','RdBu',60,'pchip'));

crange = [-1 1] * 0.15;
bins = crange(1):0.02:crange(2);

time_length = datenum(10,1,1);
yrange = [0 time_length];

xrange = [-140 -124];  % longitude range

ref_speed = 0.02;  % m/s

for i = 1:length(run)
    figure(75+i)
    clf
    set(gcf,'color','w','position',[414 126 567 826])
    
    [~,eta] = min(abs(run(i).LAT(1,:) - target_lat));
    xi = run(i).LON(:,1)>xrange(1) & run(i).LON(:,1)<xrange(2);
    
    dx = lonlat_to_dxdy(run(i).LON(1,eta), run(i).LAT(1,eta), run(i).LON(2,eta), run(i).LAT(2,eta));  % km
    x = cumsum(repmat(dx, [length(run(i).LON(xi,eta)), 1]), 'reverse');  % km
    
    ssh = run(i).ssh_raw - repmat(run(i).ssh_trend,[1,size(run(i).LON)]);  % time-mean at each point (removes large-scale structure)
    
    t_inds = find(run(i).time > run(i).time(end)-time_length);
    
    if i == 3
        t_inds = find(run(i).time > run(i).time(end)-datenum(5,1,1)-time_length);
    end
    
    time_mean = repmat(mean(ssh(:,xi,eta),1,'omitnan'), [length(t_inds),1]);
    
    time = run(i).time(t_inds)-run(i).time(t_inds(1));
    pcolor(x, time, ssh(t_inds,xi,eta) - time_mean); shading flat    

    
    colormap(jet)
    xlabel('Distance from coast (km)')
    ylabel('Time (years)')
    
    datetick('y','yy')
    
    set(gca,'fontsize',14,'ydir','normal','xdir','reverse')
    
    caxis(crange)
    ylim(yrange)
    
    colormap(cmap)

    % draw speed line
    time_sec = time*24*3600; 
    position = x(end) + (time_sec * ref_speed * 1e-3);
    line([position(1) position(end)], [time(1) time(end)], 'linewidth',4,'color',1*[1,1,1]);
    line([position(1) position(end)], [time(1) time(end)], 'linewidth',4,'color',0*[1,1,1],'linestyle','--');
    
    % draw transition line
    if i == var_res_ind
        [c.lon,c.lat] = contour_line(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour);
        transition_lon = interp1(c.lat, c.lon, target_lat);
        
        transition_x = interp1(run(i).LON(xi,eta), x, transition_lon);
        line(transition_x*[1,1], yrange, 'color',0.2*[1,1,1],'linewidth',4)
        line(transition_x*[1,1], yrange, 'color',1.0*[1,1,1],'linewidth',2)
    end    
    
    save_name = sprintf('figures/westward_propogation/hovmoeller_cal_%s_v3.png',run(i).short_name);
    set(gca,'color','none');
    export_fig(save_name, '-m3','-painters','-transparent')
end

%% colorbar
figure 
set(gcf,'color','w','position',[540 421 710 505]))
colormap(cmap)
caxis(crange)

cb = colorbar;
ylabel(cb, 'SSH anomaly (m)')
set(cb,'ydir','normal','fontsize',14)
set(gca,'visible','off')
cbarrow

save_name = 'figures/westward_propogation/hovmoeller_cal_colorbar_v3';
export_fig(save_name, '-pdf','-painters','-transparent')
