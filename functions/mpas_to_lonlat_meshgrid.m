function [LON, LAT, FIELD] = mpas_to_lonlat_meshgrid(field_to_read, mesh_fi, data_fi, lon_vec, lat_vec, t_ind)
%MPAS_TO_LONLAT_MESHGRID 
%   Detailed explanation goes here
%
%   If t_ind==0, skips time index (e.g. if reading areaCell)
%
% Kevin Rosa
% May 23, 2019


% read field from MPAS file on unstructured grid
[mpas.lon, mpas.lat] = read_mesh_file_lonlat(mesh_fi);

if t_ind>0
    mpas.field = ncread(data_fi, field_to_read, [1,t_ind], [Inf,1]);
elseif t_ind==0
    mpas.field = ncread(data_fi, field_to_read, 1, Inf);
end

% create regular grid
[LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);

% interpolate to new lon/lat matrix
dx = 4 * abs(lon_vec(2)-lon_vec(1));  % keep indices slightly larger than target region to improve interpolation near edges 
inds = mpas.lon>lon_vec(1)-dx & mpas.lon<lon_vec(end)+dx & mpas.lat>lat_vec(1)-dx & mpas.lat<lat_vec(end)+dx;
F = scatteredInterpolant(mpas.lon(inds), mpas.lat(inds), mpas.field(inds), 'linear','none'); 

FIELD = F(LON, LAT);


function [lon, lat] = read_mesh_file_lonlat(fi)
%READ_MESH_FILE_LONLAT Read lon and lat on cells and convert from radians to 180 deg coordinates

try 
    lon = rad2deg(ncread(fi, 'lonCell'));
    lat = rad2deg(ncread(fi, 'latCell'));
catch
    lon = rad2deg(ncread(fi, 'grid_center_lon'));
    lat = rad2deg(ncread(fi, 'grid_center_lat'));
end
lon(lon>180) = lon(lon>180) - 360;

end

function [LON, LAT] = make_lonlat_matrix(lon, lat)
%MAKE_LONLAT_MATRIX lon varies along LON(:,1), lat along LAT(1,:)

LON = repmat(lon(:),  [1,length(lat)]);
LAT = repmat(lat(:)', [length(lon),1]);

end

end