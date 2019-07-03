function [V] = geostrophic_vel_meridional(LON, LAT, SSH)
%GEOSTROPHIC_VEL_MERIDIONAL 
% [V] = geostrophic_vel_meridional(LON, LAT, SSH)
%
% fv = g * dn/dx  [where 'n' is free surface height]
%
% - dx taken along 1st dimension (LON(:,1) must be the lon variability).
% - LAT is necessary for calculating f
% 
% 
% Kevin Rosa
% June 11, 2019


%%
g = 9.81;  % m2/s
f = 2 * 7.2921e-5 * sind(LAT);

dlon = NaN(size(LON));
dn = NaN(size(LON));

% centered difference
dlon(2:end-1,:) = LON(3:end,:) - LON(1:end-2,:);
dx = dlon .* cosd(LAT) * (pi/180) * 6371e3;  % m

dn(2:end-1,:) = SSH(3:end,:) - SSH(1:end-2,:);

% calculate geostrophic velocity
dndx = dn ./ dx;

V = (g ./ f) .* dndx;  % m/s


end

