function [LON, LAT] = make_lonlat_matrix(lon, lat)
%MAKE_LONLAT_MATRIX lon varies along LON(:,1), lat along LAT(1,:)
% [LON, LAT] = make_lonlat_matrix(lon, lat)
%
% Kevin Rosa
% June 5, 2019

LON = repmat(lon(:),  [1,length(lat)]);
LAT = repmat(lat(:)', [length(lon),1]);

end

