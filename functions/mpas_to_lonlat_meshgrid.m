function [LON, LAT, FIELD] = mpas_to_lonlat_meshgrid(field_to_read, mesh_fi, data_fi, lon_vec, lat_vec, t_ind)
%MPAS_TO_LONLAT_MESHGRID 
%   Detailed explanation goes here
%
%   If t_ind==0, skips time index (e.g. if reading areaCell)
%   If t_ind=='all' downloads all time indices
%
% Kevin Rosa
% May 23, 2019


% read field from MPAS file on unstructured grid
[mpas.lon, mpas.lat] = read_mesh_file_lonlat(mesh_fi);

if strcmp(t_ind,'all')
    mpas.field = ncread(data_fi, field_to_read, [1,1], [Inf,Inf]);
elseif t_ind>0
    mpas.field = ncread(data_fi, field_to_read, [1,t_ind], [Inf,1]);
elseif t_ind==0
    mpas.field = ncread(data_fi, field_to_read, 1, Inf);
end

% create regular grid
[LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);

% interpolate to new lon/lat matrix
dx = 4 * abs(lon_vec(2)-lon_vec(1));  % keep indices slightly larger than target region to improve interpolation near edges 
inds = mpas.lon>lon_vec(1)-dx & mpas.lon<lon_vec(end)+dx & mpas.lat>lat_vec(1)-dx & mpas.lat<lat_vec(end)+dx;


if strcmp(t_ind,'all')  % multiple time indices
    t_length = length(mpas.field(1,:));
    FIELD = NaN([t_length, size(LON)]);
    
    for t = 1:t_length
        F = scatteredInterpolant(mpas.lon(inds), mpas.lat(inds), mpas.field(inds,t), 'linear','none'); 
        FIELD(t,:,:) = F(LON, LAT);
    end

else
    F = scatteredInterpolant(mpas.lon(inds), mpas.lat(inds), mpas.field(inds), 'linear','none');
    FIELD = F(LON, LAT);
    
end

end