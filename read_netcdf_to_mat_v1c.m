% Kevin Rosa
% July 16, 2019
%
% v1:
% - download all time indices
% v1c:
% - add arctic
% June 18, 2019

% $ cd repos/MPAS-Scripts
% $ nohup /ccs/opt/matlab-R2017b/bin/matlab -nodesktop -nosplash -r "run read_netcdf_to_mat_v1.m" &

addpath(genpath('.'))

%%
j = 1;
region(j).name = 'Pacific';
region(j).lon_range = [-143 -116];
region(j).lat_range = [30 43];
region(j).dx = 0.1;
j = j+1;
region(j).name = 'Atlantic';
region(j).lon_range = [-97 -50];
region(j).lat_range = [18 45];
region(j).dx = 0.1;
j = j+1;
region(j).name = 'NorthAtlantic';
region(j).lon_range = [-80 20];
region(j).lat_range = [45 80];
region(j).dx = 0.1;

%%
type = 'highFrequency';

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).years = 6:22;
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).years = 20:36;
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = 'EC60to30_G_case';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).years = 40:41;%56;
run(i).color = rgb('blue');

%%
FIELDS = {'kineticEnergyAtSurface','barotropicSpeed','ssh','temperatureAtSurface','temperatureAt250m','salinityAtSurface','salinityAt250m','relativeVorticityAt250m','dThreshMLD','tThreshMLD'};
t_ind = 'all';

%%
for j = 3%1:length(region)

dx = 0.1;
lon_vec = region(j).lon_range(1):region(j).dx:region(j).lon_range(2);
lat_vec = region(j).lat_range(1):region(j).dx:region(j).lat_range(2);


for i = [3,1,2]%1:length(run)
    files = {};
    for year = run(i).years
        dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.highFrequencyOutput.%04i-*',year)));
        files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
    end

    t_inds_per_year = 350/5;  % underestimate a few since it's better to be under than over
    no_years = length(files)/12;
    t_length = floor(t_inds_per_year * no_years);  % approx number of time indices -- used to initialize arrays 

    nans1d = NaN(t_length, 1);
    run(i).time = nans1d;
    
    nans3d = NaN([t_length, length(lon_vec), length(lat_vec)]);
    for F = FIELDS
        run(i).(F{1}) = nans3d;
    end

    t1 = 1;  % time index counter
    for m = 1:length(files)
        
        data_fi = files{m};
        
        new_time = mpas_time(data_fi, t_ind);
        t2 = t1 + length(new_time) - 1;
        
        run(i).time(t1:t2) = new_time;            
        
        % Calculate LON LAT and land-sea mask (only do this once)
        if t1 == 1
            [LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
            run(i).LON = LON;
            run(i).LAT = LAT;
            
            run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);
        end
        
        % read data fields
        for F = FIELDS
            [~,~,data] = mpas_to_lonlat_meshgrid(F{1}, run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);        
            run(i).(F{1})(t1:t2,:,:) = data;

        end
        
        fprintf('%s %s %.1f%s\n', region(j).name, run(i).short_name, 100*m/length(files), '%')
        t1 = t2+1;
    end
    
    % apply mask
    for F = FIELDS
        mask = repmat(permute(run(i).mask,[3,1,2]), [length(run(i).time),1,1]);
        run(i).(F{1}) = run(i).(F{1}) .* mask;
    end
    
    % save data to mat
    save_name = sprintf('%s_years%ito%i_time%s_lon%.1fto%.1f_lat_%.1fto%.1f_%s_%s.mat',...
        region(j).name, run(i).years(1), run(i).years(end), num2str(t_ind),...
        region(j).lon_range(1), region(j).lon_range(2), region(j).lat_range(1), region(j).lat_range(2),...
        type, run(i).code);
    out = run(i);
    
    save(save_name, '-struct', 'out', '-v7.3');  % -struct saves structure fields as separate variables in the mat file 
    
    % reduce memory usage a little:
    clear out
    for F = FIELDS
        run(i).(F{1}) = [];
    end
end

end
