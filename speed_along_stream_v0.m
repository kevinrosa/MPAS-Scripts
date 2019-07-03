



addpath(genpath('.'))

i = 1;
run(i).name = 'B case (fully coupled)';
run(i).short_name = 'B-case';
run(i).dir = '/scratch/kanga/runs/A_WCYC1850_ne30_oNAEC60to30cr8L60v1_anvil01/mpaso.hist.am.highFrequencyOutput/';
run(i).mesh_fi = '/scratch/kanga/runs/A_WCYC1850_ne30_oNAEC60to30cr8L60v1_anvil01/mpaso.rst.0001-01-06_00000.nc';
run(i).color = rgb('grey');

%%

files = dir(fullfile(run(i).dir, '*high*.nc'));

dx = 0.1;
lon_vec = -84:dx:-54;
lat_vec = 30:dx:42;
sz = [length(lon_vec), length(lat_vec)];

t_ind = 1;

for month = 10%:length(files)
    data_fi = fullfile(files(month).folder, files(month).name);
    [LON, LAT, ssh] = mpas_to_lonlat_meshgrid('ssh', run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
    [LON, LAT, ke_mat(:,:,1)] = mpas_to_lonlat_meshgrid('kineticEnergyAtSurface', run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
end
speed_mat = sqrt(2 * ke_mat);

run(i).lon = LON;
run(i).lat = LAT;

%%
width = 100;  % km
ssh_contour = -0.1;

[C.lon, C.lat] = streamline_coords(LON, LAT, ssh, ssh_contour);


%%
[left, right, angle] = cross_stream_transects(C.lon, C.lat, width);

%%
N = 100;
for k = 1:length(left.lon)
    flux(k) = integrate_along_transect(LON,LAT,speed_mat,left.lon(k),left.lat(k),right.lon(k),right.lat(k),N);
end

%%
dist = distance_along_stream(C.lon, C.lat);

%%
DIST = 0:20:2e3;

D.flux = interp1(dist, flux, DIST);
D.angle = interp_angle(dist, angle, DIST);

%%
figure

line(dist, flux, 'color','k')

line(DIST, D.flux, 'color','b')



%%
figure(4)
clf

pcolor(LON, LAT, speed_mat); shading flat
colormap(jet)

line(C.lon, C.lat, 'color','k')

%%
for k = 1:length(left.lon)
    pause(0.01)
    line([left.lon(k),right.lon(k)], [left.lat(k), right.lat(k)],'color','w')
end

%%

%%
figure
line(dist, flux)


