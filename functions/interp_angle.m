function angleq = interp_angle(x, angle, xq)
%INTERP_ANGLE
%   angleq = interp_angle(X, angle, Xq)
%
% Returns angleq in deg 0-360
%
% Kevin Rosa
% May 31, 2019

angle = rad2deg( unwrap(deg2rad(angle)) );  % converts to -180 180

angleq = interp1(x, angle, xq);

angleq = mod(angleq, 360);  % converts back to 0 360

end