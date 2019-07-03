% Kevin Rosa
% June 25, 2019

%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 6:22;%15:20;%6:22;
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).years = 20:36;%25:30;%20:36;
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = '20180305.GM600.T62_oECv3.eos';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).years = 20:36;%25:30;%20:36;
run(i).color = rgb('blue');


%%
bottomDepth
%%

xrange = [-85 -5];
yrange = [20 45];
dx = 0.1;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

FIELDS = {'ssh'};

for i = 1:length(run)
    files = {};
    for year = run(i).years
        dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.highFrequencyOutput.%04i-*',year)));
        files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
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
    for m = 1:length(files)
        
        data_fi = files{m};
        
        run(i).time(tt) = mpas_time(data_fi, t_ind);            

        
        % Calculate LON LAT and land-sea mask (only do this once)
        if tt == 1
            [LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
            run(i).LON = LON;
            run(i).LAT = LAT;
            
            run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);
            
            % read depths
            [~,~,run(i).bottomDepth] = mpas_to_lonlat_meshgrid('bottomDepth', run(i).mesh_fi, run(i).mesh_fi, lon_vec, lat_vec, 0);
            
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


%% SSH to geostrophic velocities
for i = 1:length(run)
    run(i).v = NaN(size(run(i).ssh));
    
    for t = 1:length(run(i).time)
        run(i).v(t,:,:) = geostrophic_vel_meridional(run(i).LON, run(i).LAT, squeeze(run(i).ssh(t,:,:)));
    end
    
    run(i).v_mean = squeeze(mean(run(i).v,1));
end

%% southward transport with m_map
% xrange = [-65 -25];
% yrange = [22 32];

figure(3)
clf
set(gcf,'position',[13 449 1629 503],'color','w')

run_inds = [3,1,2];
m = 1;
n = length(run_inds);
p = 1;

crange = [-1 1]*0.05;
dc = 0.001;
bins = 2*crange(1):dc:crange(2)*3;
% cmap = cbrewer('div','RdBu', length(bins)-1, 'pchip');
cmap = flipud(cmocean('balance',length(bins)-1));

for i = run_inds
    subplot(m,n,p)
        
    m_proj('lambert','long',xrange,'lat',yrange);
    hold on
    m_contourf(run(i).LON, run(i).LAT, run(i).v_mean, bins, 'linecolor','none'); 
%     contourf(run(i).LON, run(i).LAT, run(i).v_mean, bins, 'linecolor','none'); 
    
    m_gshhs_l('patch',0.8*[1,1,1],'edgecolor','k');
    m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')
    
%     inds = run(i).LON>=xrange(1) & run(i).LON<=xrange(2) & run(i).LAT>=yrange(1) & run(i).LAT<=yrange(2);
%     val = mean(V_mean(inds));
%     
%     % draw box
%     line([xrange(1),xrange(1),xrange(2),xrange(2),xrange(1)], [yrange(1),yrange(2),yrange(2),yrange(1),yrange(1)],'color','k','linewidth',2)
%     text(mean(xrange), mean(yrange), sprintf('mean: %.3f m/s',val),'horizontalalignment','center')
%     
    colormap(cmap)
    caxis(crange)
    
    colorbar
    title(run(i).name)
    p = p+1;
end

%%
set(gcf,'color','w')
save_name = 'figures/sverdrup/geostrophic_meridional_mean_v03.png';
saveas(gcf, save_name)


%% transport calculation
target_lat = 30;
lon_range = [-65 -20];

for i = run_inds
    [~,eta] = min(abs(run(i).LAT(1,:)-target_lat));
    xi = run(i).LON(:,eta)>=lon_range(1) & run(i).LON(:,eta)<=lon_range(2);
    
    lon = run(i).LON(xi,eta);
    lat = run(i).LAT(xi,eta);

    vel = mean(run(i).v_mean(xi,eta));
    [dx,~] = lonlat_to_dxdy(lon(1),lat(1),lon(2),lat(2));  % in km
    
    run(i).trans = sum(run(i).v_mean(xi,eta) .* run(i).bottomDepth(xi,eta) * dx*1000) * 1e-6;  % Sverdrups
    
%     fprintf('%s mean: %.3f m/s \n', run(i).short_name, mean(run(i).v_mean(xi,eta)))    
    fprintf('%s transport: %.2f Sv \n', run(i).short_name, trans)    
    
end

%% add transports to maps
p = 1;
for i = run_inds
    subplot(m,n,p)

    m_line(lon, lat, 'color','b','linewidth',2)
    txt(p) = m_text(mean(lon),mean(lat), sprintf('%.0f Sv',run(i).trans),'verticalalignment','top','horizontalalignment','center');
    
    p = p+1;
end

%% save
save_name = 'figures/sverdrup/geostrophic_meridional_mean_addtransport_v03.png';
saveas(gcf, save_name)

