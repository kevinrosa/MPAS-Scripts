function [FIELD_MEAN] = mean_along_section(field_to_read, mesh_fi, data_fi, lon_sec, lat_sec, t_ind)
%MEAN_ALONG_SECTION 
%   Detailed explanation goes here
%
%   If t_ind==0, skips time index (e.g. if reading areaCell)
%
% Kevin Rosa
% June 3, 2019

dlon = 1;  % degrees east and west of section to read mpas data for interpolation
dlat = 1;  % degrees north and south...

xrange = [min(lon_sec) - dlon, max(lon_sec) + dlon];
yrange = [max(lat_sec) - dlat, max(lat_sec) + dlat];

% make sure the lon/lat of the section are the same length (in case one is just one number) 
if length(lon_sec) == 1 
    lon_sec = lon_sec * ones(size(lat_sec));
elseif length(lat_sec) == 1 
    lat_sec = lat_sec * ones(size(lon_sec));
end

% read MPAS mesh lon/lat
[mpas.lon, mpas.lat] = read_mesh_file_lonlat(mesh_fi);

% read field (on MPAS unstructured mesh)
if t_ind>0
    mpas.field = ncread(data_fi, field_to_read, [1,t_ind], [Inf,1]);
elseif t_ind==0
    mpas.field = ncread(data_fi, field_to_read, 1, Inf);
end
% NaN out missing data
mpas.field(mpas.field<-1e4) = NaN;

% interpolate to section
inds = mpas.lon>xrange(1) & mpas.lon<xrange(end) & mpas.lat>yrange(1) & mpas.lat<yrange(end);
F = scatteredInterpolant(mpas.lon(inds), mpas.lat(inds), mpas.field(inds), 'linear','none'); 

FIELD = F(lon_sec, lat_sec);
FIELD_MEAN = mean(FIELD(:),'omitnan');

end
