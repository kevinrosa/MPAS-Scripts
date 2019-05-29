function [LON, LAT, FIELD] = mpas_to_lonlat_meshgrid(field_to_read, mesh_fi, data_fi, lon_vec, lat_vec, t_ind)
%MPAS_TO_LONLAT_MESHGRID 
%   Detailed explanation goes here
%
% Kevin Rosa
% May 23, 2019


% read field from MPAS file on unstructured grid
[mpas.lon, mpas.lat] = read_mesh_file_lonlat(mesh_fi);

mpas.field = ncread(data_fi, field_to_read, [1,t_ind], [Inf,1]);

% create regular grid
[LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);

% interpolate to new lon/lat matrix
dx = 4 * abs(lon_vec(2)-lon_vec(1));  % keep indices slightly larger than target region to improve interpolation near edges 
inds = mpas.lon>lon_vec(1)-dx & mpas.lon<lon_vec(end)+dx & mpas.lat>lat_vec(1)-dx & mpas.lat<lat_vec(end)+dx;
F = scatteredInterpolant(mpas.lon(inds), mpas.lat(inds), mpas.field(inds), 'linear','none'); 

FIELD = F(LON, LAT);


function [lon, lat] = read_mesh_file_lonlat(fi)
%READ_MESH_FILE_LONLAT Read lon and lat on cells and convert from radians to 180 deg coordinates

lon = rad2deg(ncread(fi, 'lonCell'));
lon(lon>180) = lon(lon>180) - 360;

lat = rad2deg(ncread(fi, 'latCell'));

end

function [LON, LAT] = make_lonlat_matrix(lon, lat)
%MAKE_LONLAT_MATRIX lon varies along LON(:,1), lat along LAT(1,:)

LON = repmat(lon(:),  [1,length(lat)]);
LAT = repmat(lat(:)', [length(lon),1]);

end

end