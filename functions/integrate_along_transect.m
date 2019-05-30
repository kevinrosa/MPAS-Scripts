function OUT = integrate_along_transect(LON_mat, LAT_mat, FIELD_mat, lon_0, lat_0, lon_1, lat_1, N)
%INTEGRATE_ALONG_TRANSECT
% OUT = integrate_along_transect(LON_mat, LAT_mat, FIELD_mat, lon_0, lat_0, lon_1, lat_1, N)
%
%   Numerically integrate across N evenly spaced points.
% 
% Kevin Rosa
% May 30, 2019

[tmp.dx, tmp.dy] = lonlat_to_dxdy(lon_0, lat_0, lon_1, lat_1);
width = sqrt(tmp.dx^2 + tmp.dy^2) * 1e3;  % meters
spacing = width/N;

lon_vec = linspace(lon_0, lon_1, N);
lat_vec = linspace(lat_0, lat_1, N);

field_along_transect = interp2(LON_mat', LAT_mat', FIELD_mat', lon_vec, lat_vec);

OUT = sum(field_along_transect * spacing);


end

