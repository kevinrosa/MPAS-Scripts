% Kevin Rosa
% May 16, 2019

fi = 'CUSP8/init.nc';

tmp.lon = rad2deg(ncread(fi, 'lonCell'));
tmp.lat = rad2deg(ncread(fi, 'latCell'));

tmp.lon(tmp.lon>180) = tmp.lon(tmp.lon>180) - 360;

%%
dx = 0.1;
lon = -78:dx:-61;
lat = 34:dx:43;

LON = repmat(lon(:),  [1,length(lat)]);
LAT = repmat(lat(:)', [length(lon),1]);

%%
inds = find(tmp.lon>lon(1) & tmp.lon<lon(end) & tmp.lat>lat(1) & tmp.lat<lat(end));

% total inds: 9891
% inds(end) - inds(1) = 639444  (so not much savings compared to
% downloading all 649425 indices)

%%
fi = 'CUSP8/mpaso.hist.am.highFrequencyOutput.0012-01-01_00.00.00.nc';

ssh = ncread(fi, 'pressureAdjustedSSH', [1,1], [Inf,1]);

F = scatteredInterpolant(tmp.lon, tmp.lat, ssh);

%%
M.ssh = F(LON, LAT);
M.lon = LON;
M.lat = LAT;

%%
figure

contourf(M.lon, M.lat, M.ssh)


%%
eta = 0.1;  % ssh (m)

cont = contour(M.lon, M.lat, M.ssh, eta*[1,1]);

%% create gulf stream line
C.lon = cont(1,:);
C.lat = cont(2,:);

C.lon(C.lon == eta) = NaN; 

% want to find the largest continuous contour
nan_inds = [find(isnan(C.lon)), length(C.lon)];

size_of_chunks = diff(nan_inds);
[~,I] = max(size_of_chunks);

cont_inds = (nan_inds(I)+1):(nan_inds(I+1)-1);

C.lon = C.lon(cont_inds);
C.lat = C.lat(cont_inds);

%%
line(C.lon, C.lat, 'color','k')


%% lon/lat to x/y
lon0 = median(M.lon);
lat0 = median(M.lat);

M.x = zeros(size(M.lon));
M.y = M.x;

% dx as a function of latitude
for j = 1:length(M.lat(1,:))
    dx(j) = m_lldist(M.lon(1:2,j), M.lat(1:2,j)); 
end
% dy constant
dy = m_lldist(M.lon(1,1:2), M.lat(1,1:2));

count_x = 0:(length(M.lon(:,1))-1);
count_y = 0:(length(M.lat(1,:))-1);

M.x = count_x(:) * dx(:)';  % matrix multiplication, not element-wise
M.y = repmat(dy*count_y(:)', [length(count_x),1]);

middle_ind = round(length(M.x(:,1))/2);
middle_km = M.x(middle_ind,1);

for j = 1:length(M.y(1,:))
    offset = middle_km - M.x(middle_ind,j);
    M.x(:,j) = M.x(:,j) + offset;
end

%%
figure
pcolor(M.x, M.y, M.ssh); shading flat


%% GS contour in x/y space
% simple interpolation wasn't working, will just recalculate contour
% C.x = interp1(M.lon(:), M.x(:), C.lon);
% C.y = interp1(M.lat(:), M.y(:), C.lat);

eta = -0.1;

cont = contour(M.x, M.y, M.ssh, eta*[1,1]);

% create gulf stream line
C.x = cont(1,:);
C.y = cont(2,:);

C.x(C.x == eta) = NaN; 

% want to find the largest continuous contour
nan_inds = [find(isnan(C.x)), length(C.x)];

size_of_chunks = diff(nan_inds);
[~,I] = max(size_of_chunks);

cont_inds = (nan_inds(I)+1):(nan_inds(I+1)-1);

C.x = C.x(cont_inds);
C.y = C.y(cont_inds);

%% check results
figure
pcolor(M.x, M.y, M.ssh); shading flat

line(C.x, C.y, 'color','k')


%% stream coordinates
width = 100;  % distance to each side of GS contour (km)
A = zeros(length(C.x), 2);
B = A;  s_left = A; s_right = A;
for k = 2:length(C.x)-1
    
    dx = C.x(k+1) - C.x(k-1);
    dy = C.y(k+1) - C.y(k-1);
    
    s_left(k,:) = [-dy, dx] ./ sqrt(dx^2 + dy^2);  % unit vector pointing to left of down-stream
    s_right(k,:) = -s_left(k,:);
    
    A(k,:) = s_left(k,:) * width + [C.x(k), C.y(k)];
    B(k,:) = s_right(k,:) * width + [C.x(k), C.y(k)];
    
end

%% check results
figure(16)
clf
set(gcf,'color','w')
pcolor(M.x, M.y, M.ssh); shading flat

line(C.x, C.y, 'color','k','linewidth',3)

for k = 1:5:length(A(:,1));
line([A(k,1),B(k,1)], [A(k,2),B(k,2)], 'color','w','linewidth',2)
line(A(k,1), A(k,2), 'marker','o','color','r','markerfacecolor','r')
line(B(k,1), B(k,2), 'marker','o','color','b','markerfacecolor','b')

end

colormap(jet)

set(gca,'fontsize',16)

cb = colorbar;
ylabel(cb, 'SSH (m)')
xlabel('Distance (km)')
ylabel('Distance (km)')

