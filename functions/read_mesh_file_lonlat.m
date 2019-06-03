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
