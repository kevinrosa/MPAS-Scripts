% v1: v0 had all the necessary functions, but was only written for a single
%     timestep. This script will incorporate all those pieces but apply to
%     multiple runs with multiple files each. 
%
% v2: in order to identify Gulf Stream by SSH contour, need to "de-mean"
%     the SSH. This is likely due to mass imbalance in G-cases (precip vs.
%     evap)



addpath(genpath('.'))

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
run(i).code = '20180305.GM600.T62_oECv3.eos';
run(i).dir = sprintf('/scratch/kanga/runs/%s/mpaso.hist.am.highFrequencyOutput/',run(i).code);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).years = 20:36;
run(i).color = rgb('blue');

%% get ssh mean (function of time) for each file

for i = 1:length(run)
    D = dir(sprintf('ssh_mean_lowpass*%s*mat',run(i).code));
    fi = D.name;
    run(i).S = load(fi);
end



%% Download SSH and Speed and calculate GS north wall contour
% down-stream distance (km)
DIST = 0:10:2e3;

dx = 0.1;
lon_vec = -82:dx:-54;
lat_vec = 27:dx:42;

width = 100;  % km
ssh_contour = 0.1;
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
    run(i).angle = nans;
    run(i).lon = nans;
    run(i).lat = nans;
    run(i).time = NaN(t_length, 1);
    
    run(i).speed = NaN([t_length, length(lon_vec), length(lat_vec)]);
    run(i).ssh = run(i).speed;

    tt = 1;
    for m = 1:length(files)
        
        data_fi = files{m};
        
        run(i).time(tt) = mpas_time(data_fi, t_ind);

        [LON, LAT, ssh] = mpas_to_lonlat_meshgrid('ssh', run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
        
        % Calculate land-sea mask (only do this once)
        if tt == 1
            run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);
        end
        
        ref_ssh = interp1(run(i).S.time, run(i).S.ssh_mean, run(i).time(tt));
        run(i).ssh(tt,:,:)  = ssh.*run(i).mask - ref_ssh;
        
        [LON, LAT, ke_mat]  = mpas_to_lonlat_meshgrid('kineticEnergyAtSurface', run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
        speed_mat = sqrt(2 * ke_mat);
        run(i).speed(tt,:,:) = speed_mat.*run(i).mask;

        [C.lon, C.lat] = streamline_coords(LON, LAT, squeeze(run(i).ssh(tt,:,:)), ssh_contour);
        
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
        
        tt = tt+1;
    end
end


%% Speed + SSH contour line 
t = 1;
field = 'ssh';
for i = 1:length(run)
    figure(65+i)
    clf
    set(gcf,'color','w')
    pcolor(LON, LAT, squeeze(run(i).(field)(t,:,:))); shading flat
    
    line(run(i).lon(t,:), run(i).lat(t,:), 'color','w','linewidth',2)
    
    cb = colorbar;
    ylabel(cb,'Surface speed (m/s)')
%     caxis([0 2])

    title(run(i).code,'interpreter','none')
    
%     saveas(gcf, sprintf('figures/geoid/speed_%s.png', code))
end


%% calculate standard deviation of GS transport
for i = 1:length(run)
    run(i).flux_mean = mean(run(i).flux, 1, 'omitnan');
    run(i).std = std(run(i).flux, 'omitnan');
end

%% plot mean and std transport
figure
for i = 1:length(run)
    line(run(i).dist, run(i).flux_mean, 'color',run(i).color,'linewidth',2)
    line(run(i).dist, run(i).flux_mean+run(i).std, 'color',run(i).color)
    line(run(i).dist, run(i).flux_mean-run(i).std, 'color',run(i).color)
end


%%
i = 2;
figure
line(run(i).time, run(i).flux(:,15))
%%
figure
for i = 1:2%length(run)
    hold on
    histogram(run(i).flux(:), 'facecolor',run(i).color)
end
%%






%%
i=1;
d = 50;
t = 1;
line(run(i).lon(t,d), run(i).lat(t,d), 'marker','o','color','r')


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
t = 10;
for i = 1:3
figure
% pcolor(LON, LAT, squeeze(run(i).speed(t,:,:))); shading flat
pcolor(LON, LAT, squeeze(run(i).ssh(t,:,:))); shading flat
line(run(i).lon(t,:), run(i).lat(t,:), 'color','w','linewidth',2)

end


%%
figure
line(run(i).time, run(i).ssh_mean)

%%
for i = 1:length(run)
    figure(55+i)
    clf
    set(gcf,'color','w')
    pcolor(LON, LAT, run(i).ssh); shading flat
    cb = colorbar;
    ylabel(cb,'SSH (m)')
%     caxis([-3.5 0.5])
    colormap(jet)
    if i == 1
        code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
    elseif i == 2
        code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
    elseif i == 3
        code = '20180305.GM600.T62_oECv3.eos';
    end
    title(code,'interpreter','none')
    
    saveas(gcf, sprintf('figures/geoid/ssh_%s.png', code))
end

%%
for i = 1:length(run)
    figure(65+i)
    clf
    set(gcf,'color','w')
    pcolor(LON, LAT, squeeze(run(i).speed(20,:,:))); shading flat
    cb = colorbar;
    ylabel(cb,'Surface speed (m/s)')
    caxis([0 2])
    if i == 1
        code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
    elseif i == 2
        code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
    elseif i == 3
        code = '20180305.GM600.T62_oECv3.eos';
    end
    title(code,'interpreter','none')
    
    saveas(gcf, sprintf('figures/geoid/speed_%s.png', code))
end

%%
run(1).fi = fullfile(run(1).dir, 'mpaso.hist.am.highFrequencyOutput.0010-01-01_00.00.00.nc');
run(2).fi = fullfile(run(2).dir, 'mpaso.hist.am.highFrequencyOutput.0026-01-01_00.00.00.nc');
run(3).fi = fullfile(run(3).dir, 'mpaso.hist.am.highFrequencyOutput.0026-01-01_00.00.00.nc');

for i = 1:length(run)
    [LON, LAT, temp] = mpas_to_lonlat_meshgrid('temperatureAt250m', run(i).mesh_fi, run(i).fi, lon_vec, lat_vec, t_ind);
    temp(temp<-10) = NaN;
    run(i).temp = temp;
end

%%
for i = 1:length(run)
    figure(75+i)
    pcolor(LON, LAT, run(i).temp); shading flat
    caxis([10 20])
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

%%
i = 2;
t = 40;
figure(40)
clf
pcolor(LON, LAT, squeeze(run(i).speed(t,:,:))); shading flat
line(run(i).lon(t,:), run(i).lat(t,:), 'color','w')

%%
figure
pcolor(LON, LAT, ssh); shading flat

%%
[LON, LAT, ssh2] = mpas_to_lonlat_meshgrid('pressureAdjustedSSH', run(i).mesh_fi, files{m}, lon_vec, lat_vec, t_ind);
figure
pcolor(LON, LAT, ssh2); shading flat

%%
[LON, LAT, temp] = mpas_to_lonlat_meshgrid('temperatureAt250m', run(i).mesh_fi, files{m}, lon_vec, lat_vec, t_ind);
%%
figure
temp(temp<-10) = NaN;
pcolor(LON, LAT, temp); shading flat

%%
hold on
contour(LON, LAT, temp, 15*[1,1],'color','w')





