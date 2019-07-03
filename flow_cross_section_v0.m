% Kevin Rosa
% July 3, 2019

addpath(genpath('.'))

%%
type = 'mpaso.hist.am.timeSeriesStatsMonthly';

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 2:6;%:30;
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).years = 1:4;%:19;
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = '20180305.GM600.T62_oECv3.eos';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).years = 1; %25:30;%20:36;
run(i).color = rgb('blue');

%%
MONTHS = 1:12;

xrange = [-135 -122];
yrange = [38 42];
dx = 0.08;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

FIELDS3D = {'timeMonthly_avg_velocityZonal'};%,'timeMonthly_avg_potentialDensity'};%'timeMonthly_avg_velocityMeridional'};

for i = 1:2%:length(run)
    files = {};
    for year = run(i).years
        for month = MONTHS
            dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.timeSeriesStatsMonthly.%04i-%02i-01.nc',year,month)));
            files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
        end
    end
    
    % Calculate LON LAT and land-sea mask (only do this once)
    [LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
    run(i).LON = LON;
    run(i).LAT = LAT;

    run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);

    % Determine which z-levels to read
    [~,~,run(i).H] = mpas_to_lonlat_meshgrid('bottomDepth', run(i).mesh_fi, run(i).mesh_fi, lon_vec, lat_vec, 0);

    z_level_depths = ncread(run(i).mesh_fi, 'refBottomDepth');
    depth_inds = find(z_level_depths < max(run(i).H(:)));
    depth_inds = [depth_inds; depth_inds(end)+1];

    run(i).z = z_level_depths(depth_inds);


    % Initialize matrices with NaNs:
    t_ind = 1;   % netcdf time index to read
    t_length = length(files);  % number of time indices 

    nans1d = NaN(t_length, 1);
    run(i).time = nans1d;
    
    nans3d = NaN([t_length, length(depth_inds), length(lon_vec), length(lat_vec)]);
    for F = FIELDS3D
        run(i).(F{1}) = nans3d;
    end

    tt = 1;
    for m = 1:length(files)
        
        data_fi = files{m};
        
        run(i).time(tt) = mpas_time(data_fi, t_ind, 'xtime_startMonthly');   
        
        
        % read data fields
        for F = FIELDS3D
            
            [~,~,data] = mpas_to_lonlat_meshgrid3d(F{1}, run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind, depth_inds);
            run(i).(F{1})(tt,:,:,:) = data;
        end
        
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
        tt = tt+1;
    end
    
    % apply mask
    for F = FIELDS3D
        mask = repmat(permute(run(i).mask,[4,3,1,2]), [length(run(i).time),length(depth_inds),1,1]);
        run(i).(F{1}) = run(i).(F{1}) .* mask;
    end
end

%%
i = 1;
lat_target = 40;
[~,eta] = min(abs(run(i).LAT(1,:)-lat_target));

t = 1;
field = 'timeMonthly_avg_velocityZonal';
x = squeeze(run(i).LON(:,eta));
X = repmat(x(:)', [length(run(i).z), 1]);
Y = repmat(run(i).z(:), [1, length(x)]);

D = squeeze(mean(run(i).(field)(:,:,:,eta), 1));

%%
crange = [-1 1]*0.1;
dc = 0.01;
bins = 2*crange(1):dc:crange(2)*2;
cmap = cbrewer('div','RdBu',length(crange(1):dc:crange(2))-1,'pchip');

figure
contourf(X, Y, D, bins)

set(gca,'ydir','reverse')
colormap(cmap)
colorbar
caxis(crange)

ylim([0 400])
