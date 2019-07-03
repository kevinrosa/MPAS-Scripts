function [MAG] = gradient_magnitude(LON, LAT, DATA)
% GRADIENT_MAGNITUDE
% MAG = gradient_magnitude(LON, LAT, DATA)
%
%  Centered difference calculation
%
%  See Vazquez-Cuervo et al. (2012)
% 
% 
% Kevin Rosa
% June 24, 2019


dlon = NaN(size(LON));
dlat = dlon;
grad_x = dlon;
grad_y = dlat;

% centered difference
dlon(2:end-1,:) = LON(3:end,:) - LON(1:end-2,:);
dlat(:,2:end-1) = LAT(:,3:end) - LAT(:,1:end-2);

dx = dlon .* cosd(LAT) * (pi/180) * 6371e3;  % m
dy = dlat .* (pi/180) * 6371e3; 

% gradient in each direction:
grad_x(2:end-1,:) = ( DATA(3:end,:)-DATA(1:end-2,:) ) ./ (2*dx(2:end-1,:));
grad_y(:,2:end-1) = ( DATA(:,3:end)-DATA(:,1:end-2) ) ./ (2*dy(:,2:end-1));

% magnitude
MAG = sqrt( grad_x.^2 + grad_y.^2 );


end

