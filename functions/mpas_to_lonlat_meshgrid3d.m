function [LON, LAT, FIELD] = mpas_to_lonlat_meshgrid3d(field_to_read, mesh_fi, data_fi, lon_vec, lat_vec, t_ind, depth_inds)
%MPAS_TO_LONLAT_MESHGRID3D
%
%   3d currently only handles single t_ind.
%   2D interpolates for each level in depth_inds.
%   zeros are turned into NaNs.
%
%
% Kevin Rosa
% July 3, 2019


                        
% read field from MPAS file on unstructured grid
[mpas.lon, mpas.lat] = read_mesh_file_lonlat(mesh_fi);

% read all depth inds (but will interpolate across them one at a time)
mpas.field = ncread(data_fi, field_to_read, [depth_inds(1),1,t_ind], [length(depth_inds),Inf,1]);

% NaNs
mpas.field(mpas.field==0) = NaN;

% create regular grid
[LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);

% interpolate to new lon/lat matrix
dx = 4 * abs(lon_vec(2)-lon_vec(1));  % keep indices slightly larger than target region to improve interpolation near edges 
inds = mpas.lon>lon_vec(1)-dx & mpas.lon<lon_vec(end)+dx & mpas.lat>lat_vec(1)-dx & mpas.lat<lat_vec(end)+dx;


FIELD = NaN([length(depth_inds), size(LON)]);

for k = 1:length(depth_inds)
    F = scatteredInterpolant(mpas.lon(inds), mpas.lat(inds), mpas.field(k,inds)', 'linear','none');
    FIELD(k,:,:) = F(LON, LAT);
    
end


end