% Kevin Rosa
% July 8, 2019

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
run(i).code = 'EC60to30_G_case';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).color = rgb('blue');
% i = i+1;
% run(i).name = 'Low-resolution G case';
% run(i).short_name = 'low-res';
% run(i).code = '20180305.GM600.T62_oECv3.eos';
% run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
% run(i).color = rgb('blue');

%% Settings
% download range:
xrange = [-143 -122];
yrange = [39 41];

% which .mat files to read:
target_string = 'time1_lon-143.0to-116.0_lat_30.0to43.0_highFrequency';

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
target_lat = 40;

% cmap = cbrewer('qual','Accent',60,'pchip');

% cmap = flipud(hsv(60));
% cmap(37:43,:) = [];
% cmap = flipud(cmap);

cmap = jet(60);

crange = [-0.15 0.13];
ylims = [0 datenum(17,1,1)];

for i = 1:3%length(run)
    figure(75+i)
    clf
    set(gcf,'color','w','position',[538 59 576 893])
    
    [~,eta] = min(abs(run(i).LAT(1,:) - target_lat));
    
    normalize = repmat(mean(run(i).ssh(:,:,eta),1,'omitnan'), [length(run(i).time),1]);
    imagesc(squeeze(run(i).LON(:,eta)), run(i).time-run(i).time(1), run(i).ssh(:,:,eta)-normalize)
    
    
    colormap(jet)
    xlabel('Longitude (^o East)')
    ylabel('Time (years)')
    
    datetick('y','yy')
    
    title(sprintf('SSH at %.1f^o N: %s', target_lat, run(i).name))
    set(gca,'fontsize',14)
    
    caxis(crange)
    xlim([-140 -123.5])
    ylim(ylims)
    
    colormap(cmap)
    
%     saveas(gcf,sprintf('figures/westward_propogation/imagesc_%s_normalize_00v0.png',run(i).short_name))
end


%% draw slopeline
% (used ginput)
xx = [-138.5757, -128.7349];
yy = [4619.8578, 4022.6463] + 100;

i = 2;
figure(75+i)

line(xx, yy, 'color','k','linewidth',6)
line(xx, yy, 'color','w','linewidth',3,'linestyle','--')

[dx,~] = lonlat_to_dxdy(xx(1),target_lat,xx(2),target_lat);
speed = abs(dx * 1e3 / (diff(yy)*24*3600));

fprintf('speed: %.2f km/day', abs(dx / diff(yy)))

%% draw approximate transition line
i = 1;
figure(75+i)

longitude = -131;

line(longitude*[1,1], ylims, 'color','k','linewidth',3,'linestyle','--')


%% save figures
for i = 1:3
    figure(75+i)
    save_name = sprintf('figures/westward_propogation/hovmoller_cc_%s_v1.png',run(i).short_name);
    saveas(gcf, save_name)
end

%%
H = 1000;
g = 10;
f0 = 1e-4;
R = sqrt(g*H) / f0  % deformation radius. long waves are much larger than R 

B0 = 2e-11;

c = B0 * R^2
