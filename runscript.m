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
% lon_vec = -71:dx:-61;
% lat_vec = 34:dx:43;
lon_vec = -77:dx:-59;
lat_vec = 33:dx:44;

width = 100;  % km

ssh_contour = -0.2;

%%
for i = 1:length(run)
    run(i).flux = [];
    run(i).angle = [];
    run(i).lon = [];
    run(i).lat = [];
    run(i).file_number = [];

    DIR = run(i).dir;
    D = dir(fullfile(DIR, '*high*.000*'));

    mesh_fi = fullfile(DIR, 'init.nc');

    for k = 1:length(D)
        fi = fullfile(DIR, D(k).name);

        [flux, angle, lon, lat] = cross_stream_sections(mesh_fi, fi, ssh_contour, lon_vec, lat_vec, width);

        run(i).flux = cat(1, run(i).flux, flux * 100 * 1e-6);
        run(i).angle = cat(1, run(i).angle, angle);
        run(i).lon = cat(1, run(i).lon, lon);
        run(i).lat = cat(1, run(i).lat, lat);
        
        run(i).file_number = cat(1, run(i).file_number, k*ones(size(flux)));

        fprintf('%s \n', D(k).name)
    end
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

