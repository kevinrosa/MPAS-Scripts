% Kevin Rosa
% July 9, 2019

addpath(genpath('.'))

%%
type = 'mpaso.hist.am.timeSeriesStatsMonthly';

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
% run(i).years = 20:30;
run(i).color = rgb('red');
run(i).levelp1 = 3;
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
% run(i).years = 10:19;
run(i).color = rgb('black');
run(i).levelp1 = 10;
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = 'EC60to30_G_case';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
% run(i).years = 1; %25:30;%20:36;
run(i).color = rgb('blue');

years = 'all';

%%
xrange = [-97 0];
yrange = [10 50];

dx = 0.1;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

FIELDS = {'timeMonthly_avg_velocityZonalDepthIntegrated','timeMonthly_avg_velocityMeridionalDepthIntegrated','timeMonthly_avg_waterColumnThickness'};

for i = 1:length(run)
    files = {};
%     for year = run(i).years
%         dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.timeSeriesStatsMonthly.%04i-*ubar.nc',year)));
%         files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
%     end
    dd = dir(fullfile(run(i).dir, 'mpaso.hist.am.timeSeriesStatsMonthly.*ubar.nc'));
    files = fullfile({dd(:).folder}, {dd(:).name})';

    t_ind = 1;
    t_length = length(files);  % number of time indices 

    nans1d = NaN(t_length, 1);
    run(i).time = nans1d;
    
    nans3d = NaN([t_length, length(lon_vec), length(lat_vec)]);
    for F = FIELDS
        run(i).(F{1}) = nans3d;
    end

    tt = 1;
    for m = 1:length(files)
        
        data_fi = files{m};
        
%         run(i).time(tt) = mpas_time(data_fi, t_ind, 'xtime_startMonthly');   
        STR = strsplit(data_fi, 'Monthly.');
        STR = strsplit(STR{end}, '_ubar.nc');
        run(i).time(tt) = datenum(STR{1},'yyyy-mm-dd');

        
        % Calculate LON LAT and land-sea mask (only do this once)
        if tt == 1
            [LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
            run(i).LON = LON;
            run(i).LAT = LAT;
            
            run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);
            
        end
        
        % read data fields
        for F = FIELDS
            [~,~,data] = mpas_to_lonlat_meshgrid(F{1}, run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
            run(i).(F{1})(tt,:,:) = data .* run(i).mask;
        end
        
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
        tt = tt+1;
    end
end

%% save as mat file
save(sprintf('ubar_structure_%s.mat',datestr(now,'yyyymmdd.HHMMSS')), 'run','xrange','yrange','-v7.3')


%% time-means
for i = 1:length(run)    
    run(i).ubar_mean = squeeze(mean(run(i).timeMonthly_avg_velocityZonalDepthIntegrated,1));    
    run(i).vbar_mean = squeeze(mean(run(i).timeMonthly_avg_velocityMeridionalDepthIntegrated,1));
    run(i).depth_mean = squeeze(mean(run(i).timeMonthly_avg_waterColumnThickness,1));
end

%% calculate grid cell widths
for i = 1:length(run)
    [~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
    run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;  % km
end


%%
%
%
var_res_ind = 1;
transition_contour = 15;

%% Plots of mean meridional barotropic velocity
run_inds = [3,1,2];
m = 1;
n = length(run_inds);
p = 1;

scale = 100;

crange = [-1 1]*1;
dc = 0.1;
crange = [crange(1)-dc/2, crange(2)+dc/2];
bins = crange(1):dc:crange(2);
cmap = cbrewer('div','RdBu', length(crange(1):dc:crange(2))-1, 'pchip');
% cmap = flipud(cmocean('balance',length(bins)-1));


for i = run_inds
    figure(100+i)
    clf
    set(gcf,'position',[13 472 1015 480],'color','w')
    
    m_proj('lambert','long',xrange,'lat',yrange);
    hold on
    
    data = run(i).vbar_mean * scale;
    data(data<min(bins)) = min(bins);%+dc/2;
    data(data>max(bins)) = max(bins);%-dc/2;
    
    m_contourf(run(i).LON, run(i).LAT, data, bins, 'linecolor','none'); 
    
    m_gshhs_l('patch',0.8*[1,1,1],'edgecolor','k');
    m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')
    
%     contourf(run(i).LON, run(i).LAT, run(i).vbar_mean, bins, 'linecolor','none')
    
    % add transition region contour
    if i == var_res_ind
        m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.5*[1,1,1],'linewidth',3)
    end
    
 
    colormap(cmap)
    caxis(crange)
    
    cb = colorbar;
    ylabel(cb, 'Meridional barotropic velocity (cm/s)')
%     title(run(i).name)
    p = p+1;
    
    cbarrow
    
    save_name = sprintf('figures/transport/vbar_atlantic_%s_v1.png',run(i).code);
    export_fig(save_name, '-m2')
end

%%  calculate meridional transport through section 
for target_lat = 34:-1:32
% lon_range = [-82 -8];
lon_range = [-90 0];

run_inds = [3,1,2];

fprintf('%i lat \n',target_lat)

for i = run_inds
    [~,eta] = min(abs(run(i).LAT(1,:)-target_lat));
    xi = run(i).LON(:,eta)>=lon_range(1) & run(i).LON(:,eta)<=lon_range(2);
    
    lon = run(i).LON(xi,eta);
    lat = run(i).LAT(xi,eta);

    [dx,~] = lonlat_to_dxdy(lon(1),lat(1),lon(2),lat(2));  % in kmra
    
%     run(i).trans = sum(run(i).vbar_mean(xi,eta) .* run(i).depth_mean(xi,eta) * dx*1000, 'omitnan') * 1e-6;  % Sverdrups
    vel = run(i).timeMonthly_avg_velocityMeridionalDepthIntegrated(:,xi,eta);
    depth = run(i).timeMonthly_avg_waterColumnThickness(:,xi,eta);
    
    trans = sum(vel .* depth .* dx*1e3 * 1e-6, 2, 'omitnan');  % time-varying transport (Sv) 
            
    fprintf('%s transport: %.2f +/- %.3f Sv \n', run(i).short_name, mean(trans), std(trans))   
    
    run(i).trans = trans;
    
end
fprintf('\n')

end
%% add transport section to plots and display text
for i = run_inds
    figure(100+i)
    
    xx = lon_range(1):0.1:lon_range(2);
    yy = repmat(target_lat, size(xx));
    m_line(xx, yy, 'color','k','linewidth',3)
    m_text(mean(lon_range), target_lat-0.1, sprintf('%.1f +/- %.1f Sv \n', mean(run(i).trans), std(run(i).trans)),'horizontalalignment','center','verticalalignment','top')
    
    save_name = sprintf('figures/transport/vbar_atlantic_w_transport_%s_v1.png',run(i).code);
    export_fig(save_name, '-m2')
    
end






