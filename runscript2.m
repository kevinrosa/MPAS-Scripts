% Kevin Rosa
% May 17, 2019

%%
i = 1;
run(i).name = 'CUSP8';
run(i).dir  = '/Volumes/APRICOT/CUSP8/';
i = i+1;
run(i).name = 'NA8';
run(i).dir  = '/Volumes/APRICOT/North Atlantic/';

%%

dx = 0.1;

lon_vec = -75.5:dx:-72.5;
lat_vec = 36:dx:38;

width = 100;  % km

ssh_contour = -0.2;

%%
for i = 1:length(run)
    run(i).max_speed = [];
    run(i).angle = [];
    run(i).lon = [];
    run(i).lat = [];
    run(i).file_number = [];

    DIR = run(i).dir;
    D = dir(fullfile(DIR, '*high*.000*'));

    mesh_fi = fullfile(DIR, 'init.nc');

    for k = 1:length(D)
        fi = fullfile(DIR, D(k).name);

        [max_speed, angle, lon, lat] = cross_stream_sections_max(mesh_fi, fi, ssh_contour, lon_vec, lat_vec, width);

        run(i).max_speed = cat(1, run(i).max_speed, max_speed);
        run(i).angle = cat(1, run(i).angle, angle);
        run(i).lon = cat(1, run(i).lon, lon);
        run(i).lat = cat(1, run(i).lat, lat);
        
        run(i).file_number = cat(1, run(i).file_number, k*ones(size(max_speed)));

        fprintf('%s \n', D(k).name)
    end
end

%% histogram
figure(126)

edges = [0.7:0.05:1.8];

for i = 1:length(run)
    clf
    if i == 1
        color = 'b';
    elseif i == 2
        color = 'r';
    end
    
    hold on
    histogram(run(i).max_speed, edges, 'facecolor',color);
    
    set(gca,'fontsize',16,'ytick',[],'color','none')
    xlabel('Max surface speed (m/s)')
    title(run(i).name)
    
    save_name = fullfile(sprintf('histogram_region01_%s_v00',run(i).name));
    export_fig(gcf, save_name,'-transparent','-png')
    
end


%% show box on map
x1 = lon_vec(1);
x2 = lon_vec(end);
y1 = lat_vec(1);
y2 = lat_vec(end);

xrange = [-77 -57];
yrange = [34 44];

figure(50)
clf
m_proj('lambert','long',xrange,'lat',yrange);
m_gshhs_l('patch',0.9*[1,1,1],'edgecolor','k');
m_grid_transparent('linestyle','none','linewidth',2,'tickdir','out','fontsize',14);

m_line([x1,x1,x2,x2,x1], [y1,y2,y2,y1,y1], 'color','k','linewidth',2)

set(gca,'color','none','fontsize',16)
save_name = fullfile(sprintf('map_region01_v00'));
export_fig(gcf, save_name,'-transparent','-png')


%% REGION 02

dx = 0.1;

lon_vec = -72.5:dx:-70.5;
lat_vec = 36:dx:41;

width = 100;  % km

ssh_contour = -0.2;

%%
for i = 1:length(run)
    run(i).max_speed = [];
    run(i).angle = [];
    run(i).lon = [];
    run(i).lat = [];
    run(i).file_number = [];

    DIR = run(i).dir;
    D = dir(fullfile(DIR, '*high*.000*'));

    mesh_fi = fullfile(DIR, 'init.nc');

    for k = 1:length(D)
        fi = fullfile(DIR, D(k).name);

        [max_speed, angle, lon, lat] = cross_stream_sections_max(mesh_fi, fi, ssh_contour, lon_vec, lat_vec, width);

        run(i).max_speed = cat(1, run(i).max_speed, max_speed);
        run(i).angle = cat(1, run(i).angle, angle);
        run(i).lon = cat(1, run(i).lon, lon);
        run(i).lat = cat(1, run(i).lat, lat);
        
        run(i).file_number = cat(1, run(i).file_number, k*ones(size(max_speed)));

        fprintf('%s \n', D(k).name)
    end
end

%% histogram
figure(126)

edges = [0.7:0.05:1.8];

for i = 1:length(run)
    clf
    if i == 1
        color = 'b';
    elseif i == 2
        color = 'r';
    end
    
    hold on
    histogram(run(i).max_speed, edges, 'facecolor',color);
    
    set(gca,'fontsize',16,'ytick',[],'color','none')
    xlabel('Max surface speed (m/s)')
    title(run(i).name)
    
    save_name = fullfile(sprintf('histogram_region02_%s_v00',run(i).name));
    export_fig(gcf, save_name,'-transparent','-png')
    
end

%% show box on map
x1 = lon_vec(1);
x2 = lon_vec(end);
y1 = lat_vec(1);
y2 = lat_vec(end);

xrange = [-77 -57];
yrange = [34 44];

figure(50)
clf
m_proj('lambert','long',xrange,'lat',yrange);
m_gshhs_l('patch',0.9*[1,1,1],'edgecolor','k');
m_grid_transparent('linestyle','none','linewidth',2,'tickdir','out','fontsize',14);

m_line([x1,x1,x2,x2,x1], [y1,y2,y2,y1,y1], 'color','k','linewidth',2)

set(gca,'color','none','fontsize',16)
save_name = fullfile(sprintf('map_region02_v00'));
export_fig(gcf, save_name,'-transparent','-png')



%% REGION 03

dx = 0.1;

lon_vec = -70.5:dx:-68;
lat_vec = 34.5:dx:41;

width = 100;  % km

ssh_contour = -0.2;

%%
for i = 1:length(run)
    run(i).max_speed = [];
    run(i).angle = [];
    run(i).lon = [];
    run(i).lat = [];
    run(i).file_number = [];

    DIR = run(i).dir;
    D = dir(fullfile(DIR, '*high*.000*'));

    mesh_fi = fullfile(DIR, 'init.nc');

    for k = 1:length(D)
        fi = fullfile(DIR, D(k).name);

        [max_speed, angle, lon, lat] = cross_stream_sections_max(mesh_fi, fi, ssh_contour, lon_vec, lat_vec, width);

        run(i).max_speed = cat(1, run(i).max_speed, max_speed);
        run(i).angle = cat(1, run(i).angle, angle);
        run(i).lon = cat(1, run(i).lon, lon);
        run(i).lat = cat(1, run(i).lat, lat);
        
        run(i).file_number = cat(1, run(i).file_number, k*ones(size(max_speed)));

        fprintf('%s \n', D(k).name)
    end
end

%% histogram
figure(126)

edges = [0.7:0.05:1.8];

for i = 1:length(run)
    clf
    if i == 1
        color = 'b';
    elseif i == 2
        color = 'r';
    end
    
    hold on
    histogram(run(i).max_speed, edges, 'facecolor',color);
    
    set(gca,'fontsize',16,'ytick',[],'color','none')
    xlabel('Max surface speed (m/s)')
    title(run(i).name)
    
    save_name = fullfile(sprintf('histogram_region03_%s_v00',run(i).name));
    export_fig(gcf, save_name,'-transparent','-png')
    
end

%% show box on map
x1 = lon_vec(1);
x2 = lon_vec(end);
y1 = lat_vec(1);
y2 = lat_vec(end);

xrange = [-77 -57];
yrange = [34 44];

figure(50)
clf
m_proj('lambert','long',xrange,'lat',yrange);
m_gshhs_l('patch',0.9*[1,1,1],'edgecolor','k');
m_grid_transparent('linestyle','none','linewidth',2,'tickdir','out','fontsize',14);

m_line([x1,x1,x2,x2,x1], [y1,y2,y2,y1,y1], 'color','k','linewidth',2)

set(gca,'color','none','fontsize',16)
save_name = fullfile(sprintf('map_region03_v00'));
export_fig(gcf, save_name,'-transparent','-png')



%% RESOLUTION
mesh_fi = '/Volumes/APRICOT/CUSP8/init.nc';

grd.lon = rad2deg(ncread(mesh_fi, 'lonCell'));
grd.lon(grd.lon>180) = grd.lon(grd.lon>180) - 360;

grd.lat = rad2deg(ncread(mesh_fi, 'latCell'));

grd.cell_area = ncread(mesh_fi, 'areaCell');
grd.dx = 2*sqrt(grd.cell_area/pi) / 1000;

inds = grd.lon>xrange(1) & grd.lon<xrange(2) & grd.lat>yrange(1) & grd.lat<yrange(2);

%%
F = scatteredInterpolant(grd.lon(inds), grd.lat(inds), grd.dx(inds));
LONv = (xrange(1):0.1:xrange(2));
LATv = (yrange(1):0.1:yrange(2));
LON = repmat(LONv(:), [1,length(LATv)]);
LAT = repmat(LATv(:)',[length(LONv),1]);
DX = F(LON, LAT);

%%

c = contour(LON, LAT, DX, [10,34,50]);
% contourc is a lower-level function which wouldn't plot but it can't handle matrix input 


%%
figure(50)
clf
m_proj('lambert','long',xrange,'lat',yrange);
m_gshhs_l('patch',0.9*[1,1,1],'edgecolor','k');
m_grid_transparent('linestyle','none','linewidth',2,'tickdir','out','fontsize',14);

% m_contour(grd.lon(inds), grd.lat(inds), grd.dx(inds), 34*[1,1])
hold on
% m_contour(LON, LAT, DX, 34*[1,1])
m_contour(LON, LAT, DX, [15,34],'color','k','linewidth',2)

%%
set(gca,'color','none','fontsize',16)
save_name = fullfile(sprintf('map_transition_contour_v00'));
export_fig(gcf, save_name,'-transparent','-png')















%%
clear U C
[U, C] = streamline_view(mesh_fi, fi, ssh_contour, lon_vec, lat_vec, width);

%%
figure
pcolor(U.x, U.y, U.field); shading flat
line(C.x, C.y, 'color','k', 'linewidth', 3)

for k = 1:length(C.left.x)
    line([C.left.x(k), C.right.x(k)], [C.left.y(k), C.right.y(k)])
end





%%
figure(50)
clf
for i = 1:length(run)
    figure
    for k = 1:max(run(i).file_number)
        inds = run(i).file_number == k;

        if i == 1
            color = 'b';
        elseif i == 2
            color = 'r';
        end

        line(run(i).lon(inds), run(i).lat(inds), 'color',color)

    %     pause(0.2)

    end
    
    title(run(i).name)
    
    xlim(xrange)
    ylim(yrange)
    
    set(gca,'color','none','fontsize',16)
    save_name = fullfile(run(i).dir, sprintf('histogram_%s_v00',group{1}));
    export_fig(gcf, save_name,'-transparent','-png')
end

%%
xrange = [-77 -57];
yrange = [34 44];



for i = 1:length(run)

figure(50)
clf
m_proj('lambert','long',xrange,'lat',yrange);
m_gshhs_l('patch',0.9*[1,1,1],'edgecolor','k');
% m_gshhs_l('color','k');
m_grid_transparent('linestyle','none','linewidth',2,'tickdir','out','fontsize',14);

    for k = 1:max(run(i).file_number)
        inds = run(i).file_number == k;

        if i == 1
            color = 'b';
        elseif i == 2
            color = 'r';
        end

        m_line(run(i).lon(inds), run(i).lat(inds), 'color',color)

    %     pause(0.2)

    end
    
    title(run(i).name)
    
%     xlim(xrange)
%     ylim(yrange)
    
    set(gca,'color','none','fontsize',16)
    save_name = fullfile(run(i).dir, sprintf('pathlines_%s_v00',run(i).name));
    export_fig(gcf, save_name,'-transparent','-png')
end
%% sverdrups if assume KE at surface propogates to 100 m
for i = 1:length(run)
figure(60+i)
clf

edges = 6:0.2:15;

m = 2;
n = 2;
subplot(m,n,1)
h1 = histogram(run(i).flux, edges, 'facecolor','k');
title(sprintf('%s: All', run(i).name))

subplot(m,n,3)
dtheta = 30;
inds = run(i).angle > (90-dtheta) & run(i).angle < (90+dtheta);
h2 = histogram(run(i).flux(inds), edges, 'facecolor','r');
title('Eastward')

subplot(m,n,4)
dtheta = 45;
inds = run(i).angle > (285-dtheta) & run(i).angle < (285+dtheta);
h3 = histogram(run(i).flux(inds), edges, 'facecolor','b');
title('Westward')
end
%%
for i = 1:length(run)
    for group = {'All', 'Eastward', 'Westward'}
        figure(80)
        clf

        edges = 6:0.2:15;
        
        if strcmp(group{1}, 'All')
            color = 'k';
            inds = true(size(run(i).angle));
            
        elseif strcmp(group{1}, 'Eastward')
            color = 'r';
            dtheta = 30;
            inds = run(i).angle > (90-dtheta) & run(i).angle < (90+dtheta);
            
        elseif strcmp(group{1}, 'Westward')
            color = 'b';
            dtheta = 45;
            inds = run(i).angle > (285-dtheta) & run(i).angle < (285+dtheta);
        end
        
        histogram(run(i).flux(inds), edges, 'facecolor',color);
        
        title(sprintf('%s: %s', run(i).name, group{1}))
        
        xlabel('Transport upper 100m (Sv)')
        set(gca,'fontsize',16)
        
        
        set(gca,'color','none')
        save_name = fullfile(run(i).dir, sprintf('histogram_%s_v00',group{1}));
        export_fig(gcf, save_name,'-transparent','-png')
    end
end

