function [flux, angle] = cross_stream_sections(mesh_fi, data_fi, cont)
%CROSS_STREAM_SECTIONS Integrate surface speed across sections along SSH cont
%
%
% Kevin Rosa
% May 17, 2019

% Issues:
% - the boundary rows and columns look weird. checkout imagesc(M.field').
%   this causes issues with the contour line at the edges.
% - should I control for the dx along the contour line? seems sloppy to
%   just use the line that matlab gives.
% - need a way to get back the lat/lon of the GS contour

% v2:
% - group several functions/steps into single function read_and_interp_new_field
% - change structure left(k).x to left.x(k) [same for 'right']
% - change stream-section angle so that angle indicates direction of stream
%   float (still in degrees clockwise of N
% - re-grid curvilinear x and y onto regular grid for easier interpolation
% - calculate flux through stream sections
%

%{
  1) scatteredInterp:
    lonlat (Jx1)    ->  lonlat (AxB)
    [mpas]              [meshgrid]

  2) arithmetic:
    lonlat (AxB)    ->  xy (AxB)
    [mesh]              [curvilinear]

  3) scatteredInterp:
    xy (AxB)        ->  xy (CxD)
    [curvi]             [mesh]
 
%}

dx = 0.1;
lon_vec = -71:dx:-61;
lat_vec = 34:dx:43;

t_ind = 1;

% contour that determines streamline
cont = -0.2;  % ssh (m)

width = 100;  % distance (km) on each side of contour to project to

%% Read SSH and calculate streamline and cross-stream sections
field_to_read = 'pressureAdjustedSSH';

[SSH.x, SSH.y, SSH.field] = read_and_interp_new_field(field_to_read, mesh_fi, data_fi, t_ind, lon_vec, lat_vec);

% generate contour line
[C.x, C.y] = contour_line_edgeworkaround(SSH.x, SSH.y, SSH.field, cont);

% cross-stream sections:
[left, right, angle] = cross_stream_transects(C.x, C.y, width);


%% Read surface velocity [sqrt(surfaceKE)]
field_to_read = 'kineticEnergyAtSurface';

[U.x, U.y, U.field] = read_and_interp_new_field(field_to_read, mesh_fi, data_fi, t_ind, lon_vec, lat_vec);

U.field(U.field<0) = 0;
U.field = sqrt(2 * U.field);


%% Integrate across the streamlines
N = 50;  % integration points

dx = width * 2 * 1000 / N;  % integration dx (meters)

flux = NaN(size(angle));

for k = 1:length(angle)
    
    x = linspace(left.x(k), right.x(k), N);
    y = linspace(left.y(k), right.y(k), N);
    
    u = interp2(U.x', U.y', U.field', x, y);

    flux(k) = sum(u) * dx;
end

inds = ~isnan(flux);

flux = flux(inds);
angle = angle(inds);


%% F U N C T I O N S

function [lon, lat] = read_mesh_file_lonlat(fi)
%READ_MESH_FILE_LONLAT Read lon and lat on cells and convert to 180 deg coordinates

% struct.lon = rad2deg(ncread(fi, 'lonCell'));
% struct.lat = rad2deg(ncread(fi, 'latCell'));
% 
% struct.lon(struct.lon>180) = struct.lon(struct.lon>180) - 360;

lon = rad2deg(ncread(fi, 'lonCell'));
lon(lon>180) = lon(lon>180) - 360;

lat = rad2deg(ncread(fi, 'latCell'));

end

function [LON, LAT] = make_lonlat_matrix(lon, lat)
%MAKE_LONLAT_MATRIX lon varies along LON(:,1), lat along LAT(1,:)

LON = repmat(lon(:),  [1,length(lat)]);
LAT = repmat(lat(:)', [length(lon),1]);

end


function [x, y] = lonlat_to_xy(lon, lat)
%LONLAT_TO_XY
%
% currently uses the m_lldist function

len_x = length(lon(:,1));
len_y = length(lat(1,:));

% dx as a function of latitude
for j = 1:len_y
    dx(j) = m_lldist(lon(1:2,j), lat(1:2,j)); 
end
% dy constant
dy = m_lldist(lon(1,1:2), lat(1,1:2));

% multiplied by dx (or dy). acts like a cumulative sum if repeated the dx
% (dy) values.
count_x = 0:(len_x-1);
count_y = 0:(len_y-1);

x = count_x(:) * dx(:)';  % matrix multiplication, not element-wise
y = repmat(dy*count_y(:)', [len_x,1]);

% this next part was kind of hacked together to make the section left-right
% symmetric. so the idea is to walk up each row and shift the x (lon) so
% that it's centered with all the other rows.

middle_ind = round(len_x/2);
middle_km = x(middle_ind,1);

for j = 1:len_y
    offset = middle_km - x(middle_ind,j);
    x(:,j) = x(:,j) + offset;
end

end

function [X, Y, FIELD] = curv_to_meshgrid(x, y, field)

dx = y(1,2) - y(1,1);
xrange = [min(x(:)), max(x(:))];
yrange = [min(y(:)), max(y(:))];

Nx = round(abs(diff(xrange))/dx);
Ny = round(abs(diff(yrange))/dx);

x_vec = linspace(xrange(1), xrange(2), Nx);
y_vec = linspace(yrange(1), yrange(2), Ny);

X = repmat(x_vec(:),  [1,Ny]);
Y = repmat(y_vec(:)', [Nx,1]);

F = scatteredInterpolant(x(:), y(:), field(:), 'linear', 'none');

FIELD = F(X, Y);
end



function [x, y, field] = read_and_interp_new_field(field_to_read, mesh_fi, data_fi, t_ind, lon_vec, lat_vec)

% read field from MPAS file on unstructured grid
[mpas.lon, mpas.lat] = read_mesh_file_lonlat(mesh_fi);
mpas.field = ncread(data_fi, field_to_read, [1,t_ind], [Inf,1]);

% create regular grid and interpolate field onto grid
[lon_mat, lat_mat] = make_lonlat_matrix(lon_vec, lat_vec);

% (selecting for indices in the target lat/lon range wasn't helpful for reading fewer points from netcdf
% but sub-setting the data before scatteredInterpolant solved the main bottleneck in the code: 
% that step went from 4 seconds to 0.04 seconds)
inds = mpas.lon>lon_vec(1) & mpas.lon<lon_vec(end) & mpas.lat>lat_vec(1) & mpas.lat<lat_vec(end);
F = scatteredInterpolant(mpas.lon(inds), mpas.lat(inds), mpas.field(inds)); 

curv.field = F(lon_mat, lat_mat);

% project lat/lon space onto x/y space 
[curv.x, curv.y] = lonlat_to_xy(lon_mat, lat_mat);

% grid became curvilinear when projecting onto x-y space. Now want to
% interpolate onto a uniform "meshgrid"-compliant x-y grid
[x, y, field] = curv_to_meshgrid(curv.x, curv.y, curv.field);

end


function [x, y] = contour_line(x_in, y_in, data, cont)
%CONTOUR_LINE returns a single continuous contour line

c = contour(x_in, y_in, data, cont*[1,1]);
% contourc is a lower-level function which wouldn't plot but it can't handle matrix input 

% create gulf stream line
x = c(1,:);
y = c(2,:);

x(x == cont) = NaN; 

% want to find the largest continuous contour
nan_inds = find(isnan(x));
nan_inds = [nan_inds, length(x)];  % append last point in case no NaNs or biggest section is last one. 

size_of_chunks = diff(nan_inds);
[~,I] = max(size_of_chunks);

cont_inds = (nan_inds(I)+1):(nan_inds(I+1)-1);

x = x(cont_inds);
y = y(cont_inds);
end

function [x, y] = contour_line_edgeworkaround(x_in, y_in, data, cont)
%CONTOUR_LINE returns a single continuous contour line
%
% was getting weird results from scatteredInterpolant along the edges so
% for now will remove the edges when calculating contour
%
% key piece of this function is that it removes rings
% and other features by only keeping the longest continuous contour line.
% 
% couldn't supress plotting so will just create a new figure so that it
% doesn't accidentally overwrite some important figure.
figure(999)

c = contour(x_in(2:end-1,2:end-1), y_in(2:end-1,2:end-1), data(2:end-1,2:end-1), cont*[1,1]);
% contourc is a lower-level function which wouldn't plot but it can't handle matrix input 

% create gulf stream line
x = c(1,:);
y = c(2,:);

x(x == cont) = NaN; 

% want to find the largest continuous contour
nan_inds = find(isnan(x));
nan_inds = [nan_inds, length(x)];  % append last point in case no NaNs or biggest section is last one. 

size_of_chunks = diff(nan_inds);
[~,I] = max(size_of_chunks);

cont_inds = (nan_inds(I)+1):(nan_inds(I+1)-1);

x = x(cont_inds);
y = y(cont_inds);

close 999

end

function [left, right, angle] = cross_stream_transects(x, y, width)
%CROSS_STREAM_TRANSECTS 
% Compute slope using centered difference 
% Angle is direction of streamflow in degrees clockwise of north

nans = NaN(length(x), 1);
left.x = nans;  left.y = nans;  right.x = nans;  right.y = nans; angle = nans;

for k = 2:length(x)-1
    
    dx = x(k+1) - x(k-1);
    dy = y(k+1) - y(k-1);
        
    s_left  = [-dy, dx] ./ sqrt(dx^2 + dy^2);  % unit vector pointing to left of down-stream
    s_right = -s_left;
    
    lxy = [width * s_left + [x(k), y(k)]];
    left.x(k) = lxy(1);     left.y(k) = lxy(2);
    
    rxy = width * s_right + [x(k), y(k)];
    right.x(k) = rxy(1);    right.y(k) = rxy(2);
    
    ang = rad2deg(atan2(dx, dy));
    
    if ang < 0
        ang = 360 + ang;
    end
        
    angle(k) = ang;
end
end

end
