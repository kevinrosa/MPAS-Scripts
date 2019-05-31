% v1: v0 had all the necessary functions, but was only written for a
% singlev
%     timestep. This script will incorporate all those pieces but apply to
%     multiple runs with multiple files each. 



addpath(genpath('.'))

i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).dir = '/scratch/kanga/runs/GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01/mpaso.hist.am.highFrequencyOutput/';
run(i).mesh_fi = '/scratch/kanga/runs/GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01/mpaso.rst.0001-08-01_00000.nc';
run(i).years = 10:22;
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'var-res';
run(i).dir = '/scratch/kanga/runs/20180208.GMPAS-IAF.T62_oRRS18v3.anvil/mpaso.hist.am.highFrequencyOutput/';
run(i).mesh_fi = '/scratch/kanga/runs/20180208.GMPAS-IAF.T62_oRRS18v3.anvil/ocean.oRRS18to6v3.scrip.181106.nc';
run(i).years = 20:32;
run(i).color = rgb('black');

%%
% down-stream distance (km)
DIST = 0:20:2e3;

dx = 0.1;
lon_vec = -84:dx:-54;
lat_vec = 30:dx:42;

width = 100;  % km
ssh_contour = -0.1;
N_integration = 100;  % number of cross-stream integration points

for i = 1:length(run)
    files = {};
    for year = run(i).years
        dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.highFrequencyOutput.%04i-*',year)));
        files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
    end

    t_ind = 1;
    t_length = length(files);  % number of time indices 

    run(i).dist = DIST;
    nans = NaN(t_length, length(DIST));
    run(i).flux = nans;
    run(i).lon = nans;
    run(i).lat = nans;
    run(i).time = NaN(t_length, 1);

    tt = 1;
    for m = 1:length(files)

        [LON, LAT, ssh] = mpas_to_lonlat_meshgrid('ssh', run(i).mesh_fi, files{m}, lon_vec, lat_vec, t_ind);
        [LON, LAT, ke_mat] = mpas_to_lonlat_meshgrid('kineticEnergyAtSurface', run(i).mesh_fi, files{m}, lon_vec, lat_vec, t_ind);
        speed_mat = sqrt(2 * ke_mat);

        [C.lon, C.lat] = streamline_coords(LON, LAT, ssh, ssh_contour);
        
        if ~isempty(C.lon)
            [left, right, angle] = cross_stream_transects(C.lon, C.lat, width);

            flux = NaN(size(left.lon));
            for k = 1:length(left.lon)
                flux(k) = integrate_along_transect(LON,LAT,speed_mat,left.lon(k),left.lat(k),right.lon(k),right.lat(k),N_integration);
            end

            dist = distance_along_stream(C.lon, C.lat);

            run(i).flux(tt,:) = interp1(dist, flux, DIST);
            run(i).angle(tt,:) = interp_angle(dist, angle, DIST);
            run(i).lon(tt,:) = interp1(dist, C.lon, DIST);
            run(i).lat(tt,:) = interp1(dist, C.lat, DIST);
        end

        time_str = ncread(files{m},'xtime',[1,t_ind],[Inf,1]);
        run(i).time(tt) = datenum(time_str','yyyy-mm-dd');

        tt = tt+1;
    end
end

%%
for i = 1:length(run)

    run(i).flux_mean = mean(run(i).flux, 1, 'omitnan');
    run(i).std = std(run(i).flux, 'omitnan');
end

%%
figure
for i = 1:length(run)
    line(run(i).dist, run(i).flux_mean, 'color',run(i).color)
    line(run(i).dist, run(i).flux_mean+run(i).std, 'color',run(i).color)
    line(run(i).dist, run(i).flux_mean-run(i).std, 'color',run(i).color)
end

%%
figure
line(run(i).dist, run(i).flux)

%%
figure
for i = 1:length(run)
    for m = 1:length(run(i).time)
        line(run(i).lon(m,:), run(i).lat(m,:),'color',run(i).color)
    end
    pause 
end

%%
for m = 1:length(run(i).time)
    if run(i).lat(m,20) < run(i).lat(m,2)
        line(run(i).lon(m,:), run(i).lat(m,:),'color','r')
        fprintf('uhoh\n')
%         pause
    end
end

