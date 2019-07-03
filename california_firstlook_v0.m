addpath(genpath('.'))

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

%% get ssh mean (function of time) for each file

for i = 1:length(run)
    D = dir(sprintf('ssh_mean_highfreq_lon-157*%s*mat',run(i).code));
    fi = D.name;
    run(i).S = load(fi);
end

%%
xrange = [-129 -121];
yrange = [33 42];
dx = 0.08;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

FIELDS = {'ssh','temperatureAtSurface','temperatureAt250m','relativeVorticityAt250m'};%,'dThreshMLD','tThreshMLD'};

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
        
        ref_ssh = interp1(run(i).S.time, run(i).S.ssh_mean, run(i).time(tt));
        run(i).ssh(tt,:,:)  = run(i).ssh(tt,:,:) - ref_ssh;
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
        tt = tt+1;
    end
end

%%
for fields = {'ssh','relativeVorticityAt250m','temperatureAtSurface'}
figure
clf; set(gcf,'color','w','position',[121 357 1363 594])
field = fields{1};
target_date = datenum(0,07,1);
run_inds = [3,1,2];

if strcmp('temperatureAtSurface', field)
    crange = [10 20];
elseif strcmp('relativeVorticityAt250m',field)
    crange = 1e-5 * [-1 1];
elseif strcmp('ssh',field)
    crange = 0.2 * [-1 1];
end

m = 1;
n = length(run_inds);
p = 1;  % subplot counter
for i = run_inds   
    
    
    frac_of_year = run(i).time/365.25 - floor(run(i).time/365.25);
    [~,t] = min( abs(frac_of_year - target_date/365.25));
    
    subplot(m,n,p)
    pcolor(run(i).LON, run(i).LAT, squeeze(run(i).(field)(t,:,:))); shading flat
    colorbar
    
    if strcmp('ssh',field)
        data = run(i).(field)(t,:,:);
        avg = mean(data(:),'omitnan');
    else
        avg = 0;
    end
    caxis(avg + crange)
    
    daspect([1, cosd(mean(run(i).LAT(:))), 1])
    
    title(sprintf('%s %s %s', field, run(i).short_name, datestr(run(i).time(t))))
    p = p+1;
end

% save_name = sprintf('figures/cali_001/3subplots_%s_001.png',field);
% saveas(gcf, save_name)
end

%% mean 
for fields = {'ssh'}
figure
clf; set(gcf,'color','w','position',[121 357 1363 594])
field = fields{1};
target_date = datenum(0,07,1);
run_inds = [3,1,2];

if strcmp('temperatureAtSurface', field)
    crange = [10 20];
elseif strcmp('relativeVorticityAt250m',field)
    crange = 1e-5 * [-1 1];
elseif strcmp('ssh',field)
    dc = 0.015;
    bins = -0.3:dc:0;
    crange = [bins(1) bins(end)];
end

m = 1;
n = length(run_inds);
p = 1;  % subplot counter
for i = run_inds   
    
    subplot(m,n,p)
    contourf(run(i).LON, run(i).LAT, squeeze(mean(run(i).(field), 1)), bins);
    colorbar
    
    caxis(crange)
    
    daspect([1, cosd(mean(run(i).LAT(:))), 1])
    
    title(sprintf('Mean %s %s', field, run(i).short_name))
    p = p+1;
end

save_name = sprintf('figures/cali_001/3subplots_mean_%s_cont%.3fm_001.png',field, dc);
saveas(gcf, save_name)
end

%% Standard deviation
for fields = {'ssh'}
figure
clf; set(gcf,'color','w','position',[121 357 1363 594])
field = fields{1};
target_date = datenum(0,07,1);
run_inds = [3,1,2];

if strcmp('temperatureAtSurface', field)
    crange = [10 20];
elseif strcmp('relativeVorticityAt250m',field)
    crange = 1e-5 * [-1 1];
elseif strcmp('ssh',field)
    dc = 0.0025;
    bins = 0.03:dc:0.1;
    crange = [bins(1) bins(end)];
end

m = 1;
n = length(run_inds);
p = 1;  % subplot counter
for i = run_inds   
    
    subplot(m,n,p)
    contourf(run(i).LON, run(i).LAT, squeeze(std(run(i).(field))), bins);
    colorbar
    
    caxis(crange)
    
    daspect([1, cosd(mean(run(i).LAT(:))), 1])
    
    title(sprintf('Std dev %s %s', field, run(i).short_name))
    p = p+1;
end

save_name = sprintf('figures/cali_001/3subplots_std_%s_cont%.4fm_001.png',field, dc);
saveas(gcf, save_name)
end


%% Download Grid resolution
for i = 1:length(run)
    [~,~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, lon_vec, lat_vec, 0);
    
    run(i).dx = 2 * sqrt(run(i).mask .* areaCell ./ pi) * 1e-3;  % (km)
    
end

%% Plot resolution
figure
clf; set(gcf,'color','w')%,'position',[121 357 1363 594])

run_inds = 1;%[3,1,2];

dc = 5;
bins = 5:dc:70;

m = 1;
n = length(run_inds);
p = 1;  % subplot counter
for i = run_inds   
    
    subplot(m,n,p)
    contourf(run(i).LON, run(i).LAT, run(i).dx, bins);
    colorbar
    
%     caxis(crange)
    
    daspect([1, cosd(mean(run(i).LAT(:))), 1])
    
    title(sprintf('Grid cell width (km) %s', run(i).short_name))
    p = p+1;
end

save_name = sprintf('figures/cali_001/grid_resolution_%s_cont%ikm_001.png',run(i).short_name,dc);
saveas(gcf, save_name)
    
    
    
