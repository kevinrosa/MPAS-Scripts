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
% download range:
xrange = [-97 -80];
yrange = [18 30];

% which .mat files to read:
target_string = 'time1_lon-97.0to-50.0_lat_18.0to45.0_highFrequency';

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

%% calculate streamlines
ssh_contour = 0;
run_inds = [1,2];
for i = run_inds
    for t = 1:length(run(i).time)
        [lon, lat] = streamline_coords(run(i).LON, run(i).LAT, squeeze(run(i).ssh(t,:,:)), ssh_contour);
        run(i).gs(t).lon = lon;
        run(i).gs(t).lat = lat;
       
    end
end

%% test 'er out
figure
no_t_inds = 24;
i=2;
for t = (length(run(i).time)-no_t_inds):length(run(i).time)
    line(run(i).gs(t).lon, run(i).gs(t).lat, 'color',run(i).color)
end


%%
no_t_inds = 24;
for i = run_inds

    figure(50)
    clf
    m_proj('lambert','long',xrange,'lat',yrange);
    m_gshhs_l('patch',0.8*[1,1,1],'edgecolor','k');
    m_grid_transparent('linestyle','none','linewidth',2,'tickdir','out','fontsize',14);
    
    for t = (length(run(i).time)-no_t_inds):length(run(i).time)
        m_line(run(i).gs(t).lon, run(i).gs(t).lat, 'color',run(i).color)
    end


    set(gca,'color','none','fontsize',16)
    save_name = sprintf('figures/streamlines/gs_contour%.1fm_%imonths_%s_v0',ssh_contour,no_t_inds,run(i).code);
    export_fig(gcf, save_name,'-transparent','-png')
end