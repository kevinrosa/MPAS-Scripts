function [left, right, angle] = cross_stream_transects(lon, lat, width)
%CROSS_STREAM_TRANSECTS Get coordinates on left and right of stream
%
% Compute slope using centered difference 
% Angle is direction of streamflow in degrees clockwise of north (assuming
% contour coordinates indexing increases as move downstream).
% 
%  'width' is distance on either side of stream (so 1/2 total width of
%  stream). Units are km.
%
% Kevin Rosa
% May 30, 2019

nans = NaN(length(lon), 1);
left.lon = nans;  left.lat = nans;  right.lon = nans;  right.lat = nans; angle = nans;

for k = 2:length(lon)-1
    
    [dx, dy] = lonlat_to_dxdy(lon(k-1), lat(k-1), lon(k+1), lat(k+1));
        
    s_left_hat  = [-dy, dx] ./ sqrt(dx^2 + dy^2);  % unit vector pointing to left of down-stream
    s_right_hat = -s_left_hat;
    
    % vectors with the [dx, dy] for the left and right points
    s_left  = s_left_hat * width;
    s_right = s_right_hat * width;
        
    [left.lon(k), left.lat(k)] = dxdy_to_lonlat(s_left(1), s_left(2), lon(k), lat(k));
    
    [right.lon(k), right.lat(k)] = dxdy_to_lonlat(s_right(1), s_right(2), lon(k), lat(k));
    
%     lxy = [width * s_left + [lon(k), lat(k)]];
%     left.x(k) = lxy(1);     left.y(k) = lxy(2);
%     
%     rxy = width * s_right + [lon(k), lat(k)];
%     right.x(k) = rxy(1);    right.y(k) = rxy(2);
    
    ang = rad2deg(atan2(dx, dy));
    
    if ang < 0
        ang = 360 + ang;
    end
        
    angle(k) = ang;
end
end

