function distance = distance_along_stream(lon,lat)
%DISTANCE_ALONG_STREAM in km
%   distance = distance_along_stream(lon,lat)
%
% Kevin Rosa
% May 31, 2019

% dist between points
ddist = NaN(length(lon)-1,1);

for k = 1:length(lon)-1
    [dx, dy] = lonlat_to_dxdy(lon(k),lat(k),lon(k+1),lat(k+1));
    
    ddist(k) = sqrt(dx^2 + dy^2);
    
end

distance = [0; cumsum(ddist(:))];

end

