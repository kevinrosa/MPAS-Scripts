function mask = compute_mask(mesh_fi, LON, LAT)
%COMPUTE_MASK Summary of this function goes here
%   mask = compute_mask(mesh_fi, LON, LAT)
%
%  mask: NaNs for land points, 1s for water points
%
%  Algorithm: for each point in LON LAT grid, finds nearest MPAS mesh cell.
%    If that cell is further away than the width of the cell, the point is
%    determined to be a land cell.
%
%
% Kevin Rosa
% June 4, 2019


%%
[mpas.lon, mpas.lat] = read_mesh_file_lonlat(mesh_fi);

areaCell = ncread(mesh_fi,'areaCell');
mpas.dx = 2 * sqrt(areaCell./pi) * 1e-3;  % cell diameter (km)

in_range = mpas.lon>min(LON(:)) & mpas.lon<max(LON(:)) & mpas.lat>min(LAT(:)) & mpas.lat<max(LAT(:));

% only keep cells in the lon lat range
mpas.lon = mpas.lon(in_range);
mpas.lat = mpas.lat(in_range);
mpas.dx  = mpas.dx(in_range);

%%
mask = ones(size(LON));

%%
for i = 1:length(LON(:))
    
    [dx, dy] = lonlat_to_dxdy(LON(i), LAT(i), mpas.lon, mpas.lat);

    [dist,ind] = min(sqrt(dx.^2 + dy.^2));
    
    if dist > 0.75*mpas.dx(ind)
        mask(i) = NaN;
    end
    
end


end

