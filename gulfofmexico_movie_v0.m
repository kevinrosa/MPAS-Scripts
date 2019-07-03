% Kevin Rosa
% June 26, 2019
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
% download range
xrange = [-97 -80];
yrange = [18 30];

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
        run(i).ssh_trend(t) = nanmean(data(:));
    end
    
    run(i).ssh = run(i).ssh_raw - repmat(run(i).ssh_trend,[1,size(run(i).LON)]);
end

%%
crange = [-1 1]*0.4;
dc = 0.01;
bins = 2*crange(1):dc:crange(2)*2;

cmap = flipud(cbrewer('div','Spectral',length(bins)-1,'pchip'));

%%
i = 3;

nFrames = length(run(i).time);
vidObj = VideoWriter(sprintf('mov_loopcurrent_%s_v0.avi',run(i).code));
vidObj.Quality = 100;
vidObj.FrameRate = 12;
open(vidObj)

figure
set(gcf,'color','w')

m_proj('lambert','long',xrange,'lat',yrange);
hold on
m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')

colormap(cmap)

for t = 1:length(run(i).time)
    [~,H] = m_contourf(run(i).LON, run(i).LAT, squeeze(run(i).ssh(t,:,:)), bins,'linecolor','none');
    
    caxis(crange)
    title(datestr(run(i).time(t)))
    
    writeVideo(vidObj, getframe(gcf));
    set(H,'visible','off')
end

close(vidObj);

%% 
