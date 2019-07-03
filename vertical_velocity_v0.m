% Kevin Rosa
% June 14, 2019


%%
addpath(genpath('.'))

%%
type = 'mpaso.hist.am.timeSeriesStatsMonthly';

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 2:6;
run(i).color = rgb('red');
run(i).levelp1 = 3;
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).years = 2:4;
run(i).color = rgb('black');
run(i).levelp1 = 10;
% i = i+1;
% run(i).name = 'Low-resolution G case';
% run(i).short_name = 'low-res';
% run(i).code = '20180305.GM600.T62_oECv3.eos';
% run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
% run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
% run(i).years = 25:30;%20:36;
% run(i).color = rgb('blue');

%%
MONTHS = 1:12;

% xtime_startMonthly
% timeMonthly_avg_vertVelocityTop

% refBottomDepth
% 60-layer:
%   layer 3:  30 m
%   layer 10: 100 m
%   layer 14: 140 m
% 80-layer: 
%   layer 10: 30 m
%   layer 20: 103 m
%   layer 23: 140 m
%   

xrange = [-84 -8];
yrange = [20 42];
dx = 0.1;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

FIELDS = {'timeMonthly_avg_vertVelocityTop'};

for i = 1:length(run)
    files = {};
    for year = run(i).years
        for month = MONTHS
            dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.timeSeriesStatsMonthly.%04i-%02i-01.nc',year,month)));
            files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
        end
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
    
    [mpas.lon, mpas.lat] = read_mesh_file_lonlat(run(i).mesh_fi);
    
    
    for m = 1:length(files)
        
        data_fi = files{m};
        
        run(i).time(tt) = mpas_time(data_fi, t_ind, 'xtime_startMonthly');            

        
        % Calculate LON LAT and land-sea mask (only do this once)
        if tt == 1
            [LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
            run(i).LON = LON;
            run(i).LAT = LAT;
            
            run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);
        end
        
        % read data fields
        for F = FIELDS
            
            mpas.field = squeeze(ncread(data_fi, F{1}, [run(i).levelp1,1,t_ind], [1,Inf,1]));
 
            % interpolate to new lon/lat matrix
            dx = 4 * abs(lon_vec(2)-lon_vec(1));  % keep indices slightly larger than target region to improve interpolation near edges 
            inds = mpas.lon>lon_vec(1)-dx & mpas.lon<lon_vec(end)+dx & mpas.lat>lat_vec(1)-dx & mpas.lat<lat_vec(end)+dx;
            SI = scatteredInterpolant(mpas.lon(inds), mpas.lat(inds), mpas.field(inds)', 'linear','none'); 

            FIELD = SI(LON, LAT);
            
            run(i).(F{1})(tt,:,:) = FIELD .* run(i).mask;
            
        end
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
        tt = tt+1;
    end
end


%% plotting
w_mmap = 0;

run_inds = [1,2];
m = 1;
n = length(run_inds);

version_code = 'v0';


for plot_name = {'vert_mean'}
    p = plot_name{1};
    
    
    figure
    set(gcf,'name',p,'position',[13 449 1629 503],'color','w')
    pcounter = 1;

    for i = run_inds
    subplot(m,n,pcounter)

    FIELD = 'timeMonthly_avg_vertVelocityTop';
    
    str = strsplit(p,'_');
    P.(p).operation = str{2};
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
%         if strcmp(p,'btspeed_mean') | strcmp(p,'btspeed_snap') | strcmp(p,'vort_snap')
%             contourf(run(i).LON, run(i).LAT, data, P.(p).bins,'linecolor','none')
%         else
%             contourf(run(i).LON, run(i).LAT, data, P.(p).bins,'linecolor','k')
%         end
        scale = 1e6;
        crange = [-1 1]*2;
        dc = 0.2;
        bins = -4:dc:4;
        cmap = cbrewer('div','RdBu', length(crange(1):dc:crange(2)), 'pchip');

        
        %contourf(run(i).LON, run(i).LAT, data*scale, bins,'linestyle','none') 
        
        pcolor(run(i).LON, run(i).LAT, data*scale); shading flat 
        
        caxis([-1 1]*10)
      
        
        hold on
    
        % add transition region contour
%         if i == var_res_ind
%             contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)
%         end
    end

%     caxis(crange)
    colormap(cmap)
    colorbar
    
    title(TITLE)

    pcounter = pcounter+1;
    end
    
    save_name = sprintf('figures/sverdrup/vertical_vel_mean_30m_draft0.png');%, region, p, version_code);
    saveas(gcf, save_name)
end

%%
XRANGE = [-60 -20];
YRANGE = [22 34];
xi  = run(i).LON(:,1) > XRANGE(1) & run(i).LON(:,1) < XRANGE(2);
eta = run(i).LAT(1,:) > YRANGE(1) & run(i).LAT(1,:) < YRANGE(2);

%%
for i = 1:length(run)
    data = run(i).timeMonthly_avg_vertVelocityTop(:,xi,eta);
    
    mean_vert = mean(data(:));
    fprintf('%s %.2f x10^6 m/s\n', run(i).short_name, mean_vert * 1e6);
end

% 30 m:
% var-res -0.63 x10^6 m/s
% high-res -0.68 x10^6 m/s

