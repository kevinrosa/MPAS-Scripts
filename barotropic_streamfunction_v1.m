% Kevin Rosa
% July 9, 2019

%%
addpath(genpath('.'))

%%
load('ubar_structure_20190709.152913.mat')

%%
% type = 'mpaso.hist.am.timeSeriesStatsMonthly';
% 
% i = 1;
% run(i).name = 'Coastally-refined G case';
% run(i).short_name = 'var-res';
% run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
% run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
% run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
% run(i).years = 20:30;
% run(i).color = rgb('red');
% run(i).levelp1 = 3;
% i = i+1;
% run(i).name = 'High-resolution G case';
% run(i).short_name = 'high-res';
% run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
% run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
% run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
% run(i).years = 10:19;
% run(i).color = rgb('black');
% run(i).levelp1 = 10;
% i = i+1;
% run(i).name = 'Low-resolution G case';
% run(i).short_name = 'low-res';
% run(i).code = 'EC60to30_G_case';
% run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
% run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
% run(i).years = 1; %25:30;%20:36;
% run(i).color = rgb('blue');
% 
% %%
% 
% xrange = [-97 0];
% yrange = [10 50];
% dx = 0.1;
% lon_vec = xrange(1):dx:xrange(2);
% lat_vec = yrange(1):dx:yrange(2);
% 
% FIELDS = {'timeMonthly_avg_velocityZonalDepthIntegrated','timeMonthly_avg_velocityMeridionalDepthIntegrated','timeMonthly_avg_waterColumnThickness'};
% 
% for i = 1:length(run)
%     files = {};
%     for year = run(i).years
%         dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.timeSeriesStatsMonthly.%04i-*ubar.nc',year)));
%         files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
%     end
% 
%     t_ind = 1;
%     t_length = length(files);  % number of time indices 
% 
%     nans1d = NaN(t_length, 1);
%     run(i).time = nans1d;
%     
%     nans3d = NaN([t_length, length(lon_vec), length(lat_vec)]);
%     for F = FIELDS
%         run(i).(F{1}) = nans3d;
%     end
% 
%     tt = 1;
%     for m = 1:length(files)
%         
%         data_fi = files{m};
%         
% %         run(i).time(tt) = mpas_time(data_fi, t_ind, 'xtime_startMonthly');            
% 
%         
%         % Calculate LON LAT and land-sea mask (only do this once)
%         if tt == 1
%             [LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
%             run(i).LON = LON;
%             run(i).LAT = LAT;
%             
%             run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);
%             
%         end
%         
%         % read data fields
%         for F = FIELDS
%             [~,~,data] = mpas_to_lonlat_meshgrid(F{1}, run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
%             run(i).(F{1})(tt,:,:) = data .* run(i).mask;
%         end
%         
%         
%         fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
%         tt = tt+1;
%     end
% end
% 
% 
% %%
% save(sprintf('ubar_structure_%s.mat',datestr(now,'yyyymmdd.HHMMSS')), 'run','-v7.3')

%%

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

var_res_ind = 1;
transition_contour = 15;



%% grid spacing matrices DX and DY
LON = run(i).LON;
LAT = run(i).LAT;
[dx_vec,~] = lonlat_to_dxdy(LON(1,:), LAT(1,:), LON(2,:), LAT(2,:));
DX = repmat(dx_vec(:)', [length(LON(:,1)), 1]) * 1000;
[~,dy_vec] = lonlat_to_dxdy(LON(:,1), LAT(:,1), LON(:,2), LAT(:,2));
DY = repmat(dy_vec(:), [1,length(LON(1,:))]) * 1000;

%% integrate meridional transport from west to east 
for i = 1:length(run)
    run(i).psi = NaN(size(run(i).ubar_mean));
    depth = squeeze(mean(run(i).timeMonthly_avg_waterColumnThickness, 1));

    for eta = 1:length(run(i).LON(1,:))

        run(i).psi(:,eta) = cumsum(run(i).vbar_mean(:,eta) .* depth(:,eta) .* DX(:,eta), 'omitnan') * 1e-6;

    end
end

%% plot barotropic streamfunction

bins = [-100, -30:10:-30, -25:5:30, 40, 50, 60, 100];

cmap = cmap_BlGrYeOrReVi200(length(bins)-1);

for i = 1:length(run)
    figure(300+i)
    clf
    set(gcf,'position',[13 472 1015 480],'color','w')
    
    m_proj('lambert','long',xrange,'lat',yrange);
    hold on
    
    
    data = run(i).psi;
    for b = 1:length(bins)-1
        col(b).inds = data>bins(b) & data<bins(b+1);
    end
    
    CONTS = 1:length(bins);
    VALS = mean([CONTS(1:end-1); CONTS(2:end)], 1);
    for b = 1:length(VALS)
        data(col(b).inds) = VALS(b);
    end
    
    m_contourf(run(i).LON, run(i).LAT, data .* run(i).mask, CONTS, 'color','w')
    
    m_gshhs_l('patch',0.8*[1,1,1],'edgecolor','k');
    m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')
    
    colormap(cmap)
    caxis([CONTS(1) CONTS(end)])
    
    cb = colorbar;
    ticklabels = num2cell(bins(:));
    ticklabels([1,end]) = {' '};
    set(cb,'Ticks',CONTS,'TickLabels',ticklabels)
    
    cbarrow
    
    
    

% 
%     contourf(run(i).LON, run(i).LAT, run(i).psi .* run(i).mask, bins, 'color','w'); 
% %     hold on
% %     contour(run(i).LON, run(i).LAT, run(i).psi .* run(i).mask, cont_lines, 'color','k');
%     colorbar
% 
%     caxis(crange)
%     colormap(cmap)
    
%     cbarrow

    save_name = sprintf('figures/streamfunction/bar_streamfun-v_%s_v04.png',run(i).code);
    export_fig(save_name,'-m2')
end



