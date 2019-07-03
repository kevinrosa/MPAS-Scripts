% Kevin Rosa
% June 11, 2019
%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 15:20;%6:22;
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).years = 25:30;%20:36;
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = '20180305.GM600.T62_oECv3.eos';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).years = 25:30;%20:36;
run(i).color = rgb('blue');


%%

xrange = [-80 -9];
yrange = [20 42];
dx = 0.1;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

FIELDS = {'ssh','relativeVorticityAt250m'};

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
end


%% calculate mean southward for subregion
xrange = [-65 -25];
yrange = [22 32];

figure(3)
clf
run_inds = [3,1,2];
m = 1;
n = length(run_inds);
p = 1;

cmap = cbrewer('div','RdBu', 30, 'pchip');

for i = run_inds
    subplot(m,n,p)
    
    V_mean = squeeze(mean(run(i).v,1));
    pcolor(run(i).LON, run(i).LAT, V_mean); shading flat
    
    inds = run(i).LON>=xrange(1) & run(i).LON<=xrange(2) & run(i).LAT>=yrange(1) & run(i).LAT<=yrange(2);
    val = mean(V_mean(inds));
    
    % draw box
    line([xrange(1),xrange(1),xrange(2),xrange(2),xrange(1)], [yrange(1),yrange(2),yrange(2),yrange(1),yrange(1)],'color','k','linewidth',2)
    text(mean(xrange), mean(yrange), sprintf('mean: %.3f m/s',val),'horizontalalignment','center')
    
    colormap(cmap)
    caxis(0.05*[-1 1])
    
    colorbar
    title(run(i).name)
    p = p+1;
end

%%
set(gcf,'color','w')
save_name = 'figures/sverdrup/geostrophic_meridional_mean_wbox_v02.png';
saveas(gcf, save_name)
        

%% vorticity mean

figure(4)
clf
run_inds = [3,1,2];
m = 1;
n = length(run_inds);
p = 1;

cmap = cbrewer('div','RdBu', 41, 'pchip');

for i = run_inds
    subplot(m,n,p)
    
    pcolor(run(i).LON, run(i).LAT, squeeze(mean(run(i).relativeVorticityAt250m,1))); shading flat
    
    colormap(cmap)
    caxis(1e-6*[-1 1])
    
    colorbar
    title(run(i).name)
    p = p+1;
end

set(gcf,'color','w')
save_name = 'figures/sverdrup/relvorticity250m_mean_v00.png';
saveas(gcf, save_name)


%% vorticity snapshot

figure(4)
clf
run_inds = [3,1,2];
m = 1;
n = length(run_inds);
p = 1;

cmap = cbrewer('div','RdBu', 41, 'pchip');

t = 40;

for i = run_inds
    subplot(m,n,p)
    
    pcolor(run(i).LON, run(i).LAT, squeeze(run(i).relativeVorticityAt250m(t,:,:))); shading flat
    
    colormap(cmap)
    caxis(1.5e-5*[-1 1])
    
    colorbar
    title(run(i).name)
    p = p+1;
end    

set(gcf,'color','w')
save_name = 'figures/sverdrup/relvorticity250m_snapshot_v00.png';
saveas(gcf, save_name)

