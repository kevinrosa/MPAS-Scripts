% Kevin Rosa
% July 10, 2019

%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);

%%
xrange = [-97 0];
yrange = [-5 80];

dx = 0.1;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

[LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
run(i).LON = LON;
run(i).LAT = LAT;

%% calculate grid cell widths
[~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;  % km


%% 1. empty map with lines
figure
% m_proj('Equidistant Cylindrical','long',xrange,'lat',yrange);
% 
% m_gshhs_l('patch',0.8*[1,1,1],'edgecolor','k');
% m_grid('linewi',2,'tickdir','out','fontsize',10)

m = load('functions/m_map/private/m_coasts.mat');
line(m.ncst(:,1), m.ncst(:,2),'color','k')
grid on

%%
hold on
contour(run(i).LON, run(i).LAT, run(i).widthCell, 15*[1,1],'linecolor','r')

%%
[x,y] = ginput;

%%
x = [
  -57.7397
  -56.4185
  -56.6938
  -58.8957
  -67.7033
  -71.9971
  -73.9237
  -74.4192
  -70.0153
  -68.8043
  -55.7029
  -50.3082
  -44.9135
  -44.9686
  -33.9040
  -21.5182 ];

y = [
    5.0
   10.2587
   16.6873
   20.5319
   23.6831
   24.0613
   25.2588
   29.4185
   32.6958
   36.3513
   41.6455
   41.8345
   46.4354
   49.1455
   60.5532
   65.2801 ];
   
%%
line(x, y, 'color','k','linewidth',2)

%%
for j = 1:length(x)
    fprintf('                    [\n')
    fprintf('                        %.3f,\n',x(j))
    fprintf('                        %.3f\n',y(j))
    fprintf('                    ],\n')
end


%% adding extra points to straight line transects
N = 60;  % number of points for each transect

% x1, y1, x2, y2
XY = [
    -70, 60, 6, 60
    -54, 48, -3, 48
    -62, 45, -1, 45
    -74, 40, -9, 40
    -76, 35, -6.5, 35
    -81, 30, -10, 30
    -80, 26, -15, 26
    -83, 15, -17, 15];

% XY = [-80.05, 26.5, -14, 26.5];

XY = [-48-15/60, 60.5, -55.5, 53+40/60];  % https://www.whoi.edu/page.do?pid=24635&print=this
N = 8;

for t = 1:length(XY(:,1))
    lat_ref = XY(t,2);
    
    lon = linspace(XY(t,1), XY(t,3), N);
    lat = linspace(XY(t,2), XY(t,4), N);
    
    indent = '            ';
    fprintf('%s{\n', indent)
    fprintf('%s"type": "Feature",\n', indent)
    fprintf('%s"properties": {\n', indent)
    fprintf('%s    "name": "Atlantic Crossing at %iN",\n', indent, lat_ref)
    fprintf('%s    "tags": "standard_transport_sections",\n', indent)
    fprintf('%s    "object": "transect",\n', indent)
    fprintf('%s    "component": "ocean",\n', indent)
    fprintf('%s    "author": "Kevin Rosa",\n', indent)
    fprintf('%s    "depth": "6000.0",\n', indent)
    fprintf('%s    "note": " ",\n', indent)
    fprintf('%s    "history": "Created July 15, 2019;"\n', indent)
    fprintf('%s},\n', indent)
    fprintf('%s"geometry": {\n', indent)
    fprintf('%s    "type": "LineString",\n', indent)
    fprintf('%s    "coordinates": [\n', indent)
    
    for j = 1:length(lon)
        fprintf('                    [\n')
        fprintf('                        %.3f,\n', lon(j))
        fprintf('                        %.3f\n', lat(j))
        if j ~= length(lon)
            fprintf('                    ],\n')
        elseif j == length(lon)
            fprintf('                    ]\n')
        end
    end
    
    fprintf('                ]\n')
    fprintf('            }\n') 
    fprintf('        },\n')
    
end

    


