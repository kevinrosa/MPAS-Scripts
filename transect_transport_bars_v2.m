% Kevin Rosa
% July 29, 2019
%%
addpath(genpath('.'))

%%
type = 'transects';
mask_version = 'northatlantic_06';
i = 1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = 'EC60to30_G_case';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).mask_fi = fullfile(run(i).dir, sprintf('*%s*.nc',mask_version));
run(i).color = rgb('blue');
i = i+1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).mask_fi = fullfile(run(i).dir, sprintf('*%s*.nc',mask_version));
run(i).color = rgb('black');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).mask_fi = fullfile(run(i).dir, sprintf('*%s*.nc',mask_version));
run(i).color = rgb('red');


for i = 1:length(run)
    D = dir(run(i).mask_fi);
    run(i).mask_fi = fullfile(D(1).folder, D(1).name);
end

%%
transect(1).name = 'AtlanticCrossingat26.5N';

%% Read velocities for each transect
for i = 1:length(run)
for j = 1:length(transect)
    
    fi = fullfile(run(i).dir, sprintf('transect.%s.nc',transect(j).name));

    run(i).T(j).fi = fi;
    for fields = {'time','dx','dz','velocity'}
        run(i).T(j).(fields{1}) = ncread(run(i).T(j).fi, fields{1});
    end

    run(i).T(j).mask = ones(size(run(i).T(j).velocity(:,:,1)));
    run(i).T(j).mask(run(i).T(j).velocity(:,:,1)==0) = NaN;    
    
end
end

%% repmat x and z
for i = 1:length(run)
for j = 1:length(transect)
    run(i).T(j).DX = repmat(run(i).T(j).dx',[length(run(i).T(j).dz),1]);
    run(i).T(j).DZ = repmat(run(i).T(j).dz, [1,length(run(i).T(j).dx)]);
    
    run(i).T(j).z = cumsum(run(i).T(j).dz);
    run(i).T(j).Z = repmat(run(i).T(j).z, [1,length(run(i).T(j).dx)]);
end
end

%% Determine lon/lat coordinates for each transect
for i = 1:length(run)
    
    transect_list = ncread(run(i).mask_fi, 'transectNames');
    
    verticesOnEdge = ncread(run(i).mesh_fi, 'verticesOnEdge');

    lonVertex = rad2deg(ncread(run(i).mesh_fi, 'lonVertex'));
    lonVertex(lonVertex>180) = lonVertex(lonVertex>180) - 360;
    latVertex = rad2deg(ncread(run(i).mesh_fi, 'latVertex'));
    
    for j = 1:length(transect)
        k = [];  % transect index to read
        for t = 1:length(transect_list(1,:))
            str_to_find = transect_list(:,t)';
            str_to_find(isspace(str_to_find)) = [];
            if ~isempty(strfind(str_to_find, transect(j).name))
                k = t;
            end
        end
        
        run(i).T(j).transect_index = k;
        
        edge_inds = ncread(run(i).mask_fi, 'transectEdgeGlobalIDs', [1,k], [Inf,1]);
        edge_inds(edge_inds==0) = [];
        vert_inds = verticesOnEdge(:, edge_inds);
        
        run(i).T(j).lon = lonVertex(vert_inds(2,:));
        run(i).T(j).lat = latVertex(vert_inds(2,:));
    
    end
end

%%
j = 1;
for i = 1:length(run)
    figure
    pcolor(repmat(run(i).T(j).lon',[length(run(i).T(j).dz),1]), run(i).T(j).Z, -mean(run(i).T(j).velocity,3)); shading flat
    colorbar
    caxis([-0.1 0.1])
    set(gca,'ydir','reverse')
    
end

%% Compute transports for different currents flowing through transect
j = 1;  % 26.5N

llon = -73;

% swbc
fprintf('SWBC:\n')
for i = 1:length(run)
    inds_x = run(i).T(j).lon < llon;
    inds_z = run(i).T(j).z < 1000;
    
    inds_t = true(size(run(i).T(j).time));
    inds_t = find(run(i).T(j).time > 5);
    
    trans_zxt = -1 * run(i).T(j).velocity(inds_z,inds_x,inds_t) .* repmat(run(i).T(j).DX(inds_z,inds_x), [1,1,length(inds_t)]) .* repmat(run(i).T(j).DZ(inds_z,inds_x), [1,1,length(inds_t)]);
    
    run(i).T(j).trans_swbc = squeeze(sum(sum(trans_zxt,1,'omitnan'),2,'omitnan')) * 1e-6;
    
    fprintf('%s: %.2f +/- %.2f Sv\n', run(i).short_name, mean(run(i).T(j).trans_swbc), std(run(i).T(j).trans_swbc))

end
fprintf('\n')

% dwbc
fprintf('DWBC:\n')
for i = 1:length(run)
    inds_x = run(i).T(j).lon < llon;
    inds_z = run(i).T(j).z >= 1000;
    
    inds_t = true(size(run(i).T(j).time));
    inds_t = find(run(i).T(j).time > 5);
    
    trans_zxt = -1 * run(i).T(j).velocity(inds_z,inds_x,inds_t) .* repmat(run(i).T(j).DX(inds_z,inds_x), [1,1,length(inds_t)]) .* repmat(run(i).T(j).DZ(inds_z,inds_x), [1,1,length(inds_t)]);
    
    run(i).T(j).trans_dwbc = squeeze(sum(sum(trans_zxt,1,'omitnan'),2,'omitnan')) * 1e-6;
    
    fprintf('%s: %.2f +/- %.2f Sv\n', run(i).short_name, mean(run(i).T(j).trans_dwbc), std(run(i).T(j).trans_dwbc))

end
fprintf('\n')

% rest
fprintf('rest:\n')
for i = 1:length(run)
    inds_x = run(i).T(j).lon >= llon;
    inds_z = run(i).T(j).z > -100;
    
    inds_t = true(size(run(i).T(j).time));
    inds_t = find(run(i).T(j).time > 5);
    
    trans_zxt = -1 * run(i).T(j).velocity(inds_z,inds_x,inds_t) .* repmat(run(i).T(j).DX(inds_z,inds_x), [1,1,length(inds_t)]) .* repmat(run(i).T(j).DZ(inds_z,inds_x), [1,1,length(inds_t)]);
    
    run(i).T(j).trans_rest = squeeze(sum(sum(trans_zxt,1,'omitnan'),2,'omitnan')) * 1e-6;
    
    fprintf('%s: %.2f +/- %.2f Sv\n', run(i).short_name, mean(run(i).T(j).trans_rest), std(run(i).T(j).trans_rest))

end
fprintf('\n')

%% bar plots

for F = {'trans_swbc','trans_dwbc','trans_rest'}
    figure
    set(gcf,'color','w','position',[318 14 240 938])
    field = F{1};
    hold on
    
    for i = 1:length(run)
        data = run(i).T(j).(field);
        h(i) = bar(i, mean(data), 'facecolor', run(i).color, 'barwidth',0.9);
        e(i) = errorbar(i, mean(data), std(data),'color',0.7*[1,1,1],'linewidth',2);
        e(i) = errorbar(i, mean(data), std(data),'color',rgb('goldenrod'),'linewidth',2);
    end
    ylim([-40 60])
    
    save_name = sprintf('figures/transport/bars_26.5N_%s_%.1fE_v1', field, llon);
    set(gca,'color','none','visible','off');
    export_fig(save_name, '-pdf', '-painters','-transparent','-nocrop')
end


%% read bottom depth
i = 1;
xrange = [-88 -8];
yrange = [14 44];
lon_vec = xrange(1):0.1:xrange(2);
lat_vec = yrange(1):0.1:yrange(2);

[LON,LAT,depth] = mpas_to_lonlat_meshgrid('bottomDepth', run(i).mesh_fi, run(i).mesh_fi, lon_vec, lat_vec, 0); 

%% map view
figure
set(gcf,'color','w','position',[40 297 1166 655])

bins = 0:500:6000;
m_proj('lambert','long',xrange,'lat',yrange);
hold on
m_contourf(LON, LAT, depth, bins, 'linestyle','none')
% m_contour(LON, LAT, depth, 1000:1000:6000, 'linecolor','k')
m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linestyle','none','box','fancy','tickdir','none','xticklabels','','yticklabels','')

% show transect line
m_line(linspace(run(i).T(j).lon(1),run(i).T(j).lon(end),300), linspace(run(i).T(j).lat(1),run(i).T(j).lat(end),300), 'color','k','linewidth',3)

cmap = cbrewer('seq','Greys',length(bins)*4,'pchip');
cmap = cmap(1:length(bins)-1, :);
colormap(cmap)

save_name = sprintf('figures/transport/map_bathy_line_26.5N_v2');
set(gca,'color','none');
export_fig(save_name, '-pdf', '-painters','-transparent')

%% cross-section
figure
set(gcf,'color','w','position',[560 516 799 410])

lat = 26.5;
[~,eta] = min(abs(LAT(1,:)-lat));

line(LON(:,eta), depth(:,eta), 'color','k','linewidth',2)

% line(-74*[1,1], [0 6000], 'color','b','linewidth',3,'linestyle','--')
% line([-90, -74], 1000*[1,1], 'color','b','linewidth',3,'linestyle','--')

xlim([-82 -11])
set(gca,'ydir','reverse')
xlabel('Longitude (^o East)')
ylabel('Bottom depth (m)')
set(gca,'fontsize',18)


save_name = sprintf('figures/transport/cross-sec_bathy_line_%.1fN_v2fs18',lat);
% set(gca,'color','none');
export_fig(save_name, '-pdf', '-painters')%,'-transparent')

%%
% save_name = sprintf('figures/transport/cross-sec_bathy_line_30N_background_v0.png');
% set(gca,'color','w');
% export_fig(save_name, '-m3','-painters','-transparent')








