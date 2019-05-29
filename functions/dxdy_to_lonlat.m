function [lon_1, lat_1] = dxdy_to_lonlat(dx, dy, lon_0, lat_0)
%LONLAT_TO_DXDY Summary of this function goes here
%   Given: dx and dy in km and lon_0 and lat_0 in degrees
%   Returns: lon_1 and lat_1 in degrees
%
% Kevin Rosa
% May 23, 2019

R = 6371;  % radius of Earth in km
lat_ref = lat_0;

A = (pi/180) * R * cosd(lat_ref);  % dist (km) per degree longitude
B = (pi/180) * R;  % dist (km) per degree latitude

dlon = dx / A;
dlat = dy / B;

lon_1 = lon_0 + dlon;
lat_1 = lat_0 + dlat;

end

