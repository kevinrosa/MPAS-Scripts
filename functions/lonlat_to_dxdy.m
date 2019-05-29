function [dx, dy] = lonlat_to_dxdy(lon_0, lat_0, lon_1, lat_1)
%LONLAT_TO_DXDY Summary of this function goes here
%   Given: lon and lat in degrees
%   Returns: dx and dy in km
%
% Kevin Rosa
% May 23, 2019

R = 6371;  % radius of Earth in km
lat_ref = mean([lat_0, lat_1]);

A = (pi/180) * R * cosd(lat_ref);  % dist (km) per degree longitude
B = (pi/180) * R;  % dist (km) per degree latitude

dx = (lon_1-lon_0) * A;
dy = (lat_1-lat_0) * B;

end

