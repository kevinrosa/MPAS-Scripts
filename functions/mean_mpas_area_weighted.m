function average = mean_mpas_area_weighted(field_to_read, mesh_fi, data_fi, lon_range, lat_range, t_ind)
%MEAN_MPAS_AREA_WEIGHTED
% average = mean_mpas_area_weighted(field_to_read, mesh_fi, data_fi, lon_range, lat_range, t_ind)
%
% Kevin Rosa
% June 4, 2019


%%
[mpas.lon, mpas.lat] = read_mesh_file_lonlat(mesh_fi);

mpas.areaCell = ncread(mesh_fi,'areaCell');

in_range = mpas.lon>lon_range(1) & mpas.lon<lon_range(end) & mpas.lat>lat_range(1) & mpas.lat<lat_range(end);

% only keep cells in the lon lat range
mpas.lon = mpas.lon(in_range);
mpas.lat = mpas.lat(in_range);
mpas.areaCell  = mpas.areaCell(in_range);

mpas.ssh = ncread(data_fi, field_to_read, [1,t_ind], [Inf,1]);
mpas.ssh = mpas.ssh(in_range,:);

%%
integral = sum(mpas.areaCell .* mpas.ssh);
average = integral / sum(mpas.areaCell);


end

