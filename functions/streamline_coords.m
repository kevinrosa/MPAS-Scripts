function [contour_lon, contour_lat] = streamline_coords(LON, LAT, DATA, CONTOUR)
%STREAMLINE_COORDS Returns longest continuous contour line
% [contour_lon, contour_lat] = streamline_coords(LON, LAT, DATA, CONTOUR)
%  
%   LON, LAT, and DATA must be a transposed meshgrid (DATA' is a meshgrid).
%   
%   Will likely change other functions soon that everything is a normal
%   Matlab meshgrid. 
%   
%   This function makes use of the low-level contourc() function. It's 
%   faster than contour() and suppresses plotting but it's less flexible
%   (can't handle curvilinear grids).
%
% Kevin Rosa
% May 23, 2019

c = contourc(LON(:,1), LAT(1,:), DATA', CONTOUR*[1,1]);

% create gulf stream line
x = c(1,:);
y = c(2,:);

x(x == CONTOUR) = NaN; 

% want to find the largest continuous contour
nan_inds = find(isnan(x));
nan_inds = [nan_inds, length(x)];  % append last point in case no NaNs or biggest section is last one. 

size_of_chunks = diff(nan_inds);
[~,I] = max(size_of_chunks);

cont_inds = (nan_inds(I)+1):(nan_inds(I+1)-1);

contour_lon = x(cont_inds);
contour_lat = y(cont_inds);

end

