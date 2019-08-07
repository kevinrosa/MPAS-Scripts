% Kevin Rosa
% July 12, 2019
%%
addpath(genpath('.'))

%%
type = 'transects';
mask_version = 'northatlantic_05';
i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).mask_fi = fullfile(run(i).dir, sprintf('*%s*.nc',mask_version));
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).mask_fi = fullfile(run(i).dir, sprintf('*%s*.nc',mask_version));
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = 'EC60to30_G_case';
run(i).dir = sprintf('/scratch/kanga/runs/%s/%s/',run(i).code,type);
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).mask_fi = fullfile(run(i).dir, sprintf('*%s*.nc',mask_version));
run(i).color = rgb('blue');

for i = 1:length(run)
    D = dir(run(i).mask_fi);
    run(i).mask_fi = fullfile(D(1).folder, D(1).name);
end

%%
transect(1).name = 'AtlanticCrossingat45N';
transect(2).name = 'AtlanticCrossingat30N';
transect(3).name = 'AtlanticCrossingat26N';
transect(4).name = 'AtlanticCrossingat15N';
transect(5).name = 'Florida-Cuba';

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

%%
for i = 1:length(run)
for j = 1:length(transect)
    run(i).T(j).DX = repmat(run(i).T(j).dx',[length(run(i).T(j).dz),1]);
    run(i).T(j).DZ = repmat(run(i).T(j).dz, [1,length(run(i).T(j).dx)]);
    
    run(i).T(j).z = cumsum(run(i).T(j).dz);
    run(i).T(j).Z = repmat(run(i).T(j).z, [1,length(run(i).T(j).dx)]);
end
end

%% Determine coordinates for each transect
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
j = 2;
for i = 1:length(run)
    figure
    
%     X = repmat(cumsum(run(i).T(j).dx)', [length(run(i).T(j).dz),1]);
    X = repmat(run(i).T(j).lon', [length(run(i).T(j).dz),1]);
    Z = repmat(cumsum(run(i).T(j).dz), [1,length(run(i).T(j).dx)]);
    
    pcolor(X, Z, mean(run(i).T(j).velocity,3) .* run(i).T(j).mask); shading flat

    set(gca,'ydir','reverse')

    crange = [-1 1] * 0.03;
    caxis(crange)
    colormap(cbrewer('div','RdBu',40,'pchip'))
    
end

%%
figure
m = load('functions/m_map/private/m_coasts.mat');
line(m.ncst(:,1), m.ncst(:,2),'color','k')
grid on

for i = 1:length(run)
    line(run(i).T(j).lon, run(i).T(j).lat, 'color',run(i).color)

end

%%
j = 2;
for i = 1:length(run)
    figure

    X = repmat(run(i).T(j).lon', [length(run(i).T(j).dz),1]);
    Z = repmat(cumsum(run(i).T(j).dz), [1,length(run(i).T(j).dx)]);

    pcolor(X, Z, mean(run(i).T(j).velocity,3) .* run(i).T(j).mask); shading flat

    set(gca,'ydir','reverse')

    crange = [-1 1] * 0.03;
    caxis(crange)
    colormap(cbrewer('div','RdBu',40,'pchip'))
end


%% computing transports: dwbc
j = 2;  % 30N

for i = 1:length(run)
    dwbc_x = run(i).T(j).lon > -76.7 & run(i).T(j).lon < -75;
    
    dwbc_x = run(i).T(j).lon > -77 & run(i).T(j).lon < -74;
       
    dwbc_z = run(i).T(j).z > 0;%700;
    
    vel = mean(run(i).T(j).velocity,3);
    transport_xz = vel(dwbc_z, dwbc_x) .* run(i).T(j).DX(dwbc_z, dwbc_x) .* run(i).T(j).DZ(dwbc_z, dwbc_x);
    
    fprintf('%s: %.2f Sv\n', run(i).short_name, sum(transport_xz(:) * 1e-6,'omitnan'))

end
fprintf('\n')

%% computing transports: sverdrup
j = 2;  % 30N

for i = 1:length(run)
    inds_x = run(i).T(j).lon > -74;
    inds_z = true(size(run(i).T(j).z));
    
    vel = mean(run(i).T(j).velocity,3);
    transport_xz = vel(inds_z, inds_x) .* run(i).T(j).DX(inds_z, inds_x) .* run(i).T(j).DZ(inds_z, inds_x);
    
    fprintf('%s: %.2f Sv\n', run(i).short_name, sum(transport_xz(:) * 1e-6,'omitnan'))

end
fprintf('\n')


%% computing transports: full transect
j = 2;  % 30N

for i = 1:length(run)
    inds_x = true(size(run(i).T(j).lon));
    inds_z = true(size(run(i).T(j).z));
    
    vel = mean(run(i).T(j).velocity,3);
    transport_xz = vel(inds_z, inds_x) .* run(i).T(j).DX(inds_z, inds_x) .* run(i).T(j).DZ(inds_z, inds_x);
    
    fprintf('%s: %.2f Sv\n', run(i).short_name, sum(transport_xz(:) * 1e-6,'omitnan'))

end
fprintf('\n')

%% computing transports: gulf stream
j = 2;  % 30N

for i = 1:length(run)
    inds_x = run(i).T(j).lon < -78.5;
    inds_z = true(size(run(i).T(j).z));
    
    vel = mean(run(i).T(j).velocity,3);
    transport_xz = vel(inds_z, inds_x) .* run(i).T(j).DX(inds_z, inds_x) .* run(i).T(j).DZ(inds_z, inds_x);
    
    fprintf('%s: %.2f Sv\n', run(i).short_name, sum(transport_xz(:) * 1e-6,'omitnan'))

end
fprintf('\n')

%% computing transports: swbc
j = 2;  % 30N

for i = 1:length(run)
    inds_x = run(i).T(j).lon < -74;
    inds_z = run(i).T(j).z < 1000;
    
    vel = mean(run(i).T(j).velocity,3);
    transport_xz = vel(inds_z, inds_x) .* run(i).T(j).DX(inds_z, inds_x) .* run(i).T(j).DZ(inds_z, inds_x);
    
    fprintf('%s: %.2f Sv\n', run(i).short_name, sum(transport_xz(:) * 1e-6,'omitnan'))

end
fprintf('\n')

%% computing transports: dwbc 2
j = 2;  % 30N

for i = 1:length(run)
    inds_x = run(i).T(j).lon < -74;
    inds_z = run(i).T(j).z >= 1000;
    
    vel = mean(run(i).T(j).velocity,3);
    transport_xz = vel(inds_z, inds_x) .* run(i).T(j).DX(inds_z, inds_x) .* run(i).T(j).DZ(inds_z, inds_x);
    
    fprintf('%s: %.2f Sv\n', run(i).short_name, sum(transport_xz(:) * 1e-6,'omitnan'))

end

%% computing transports: florida-cuba
j = 5; % florida-cuba
for i = 1:length(run)
    inds_x = true(size(run(i).T(j).lon));
    inds_z = true(size(run(i).T(j).z));
    
    vel = mean(run(i).T(j).velocity,3);
    transport_xz = vel(inds_z, inds_x) .* run(i).T(j).DX(inds_z, inds_x) .* run(i).T(j).DZ(inds_z, inds_x);
    
    fprintf('%s: %.2f Sv\n', run(i).short_name, sum(transport_xz(:) * 1e-6,'omitnan'))

end
fprintf('\n')

%%
%
%
%% computing transports: sverdrup
j = 3;  % 26N

for i = 1:length(run)
    inds_x = run(i).T(j).lon > -73;
    inds_z = true(size(run(i).T(j).z));
    
    vel = mean(run(i).T(j).velocity,3);
    transport_xz = vel(inds_z, inds_x) .* run(i).T(j).DX(inds_z, inds_x) .* run(i).T(j).DZ(inds_z, inds_x);
    
    fprintf('%s: %.2f Sv\n', run(i).short_name, sum(transport_xz(:) * 1e-6,'omitnan'))

end
fprintf('\n')
%% computing transports: swbc
j = 3;  % 26N

for i = 1:length(run)
    inds_x = run(i).T(j).lon < -73;
    inds_z = run(i).T(j).z < 1000;
    
    vel = mean(run(i).T(j).velocity,3);
    transport_xz = vel(inds_z, inds_x) .* run(i).T(j).DX(inds_z, inds_x) .* run(i).T(j).DZ(inds_z, inds_x);
    
    fprintf('%s: %.2f Sv\n', run(i).short_name, sum(transport_xz(:) * 1e-6,'omitnan'))

end
fprintf('\n')

%% computing transports: dwbc
j = 3;  % 26N

for i = 1:length(run)
    inds_x = run(i).T(j).lon < -73;
    inds_z = run(i).T(j).z >= 1000;
    
    vel = mean(run(i).T(j).velocity,3);
    transport_xz = vel(inds_z, inds_x) .* run(i).T(j).DX(inds_z, inds_x) .* run(i).T(j).DZ(inds_z, inds_x);
    
    fprintf('%s: %.2f Sv\n', run(i).short_name, sum(transport_xz(:) * 1e-6,'omitnan'))

end

%%

%%
fi = '/scratch/kanga/runs/EC60to30_G_case/transects/oEC60to30v3_sections_northatlantic_03.nc';
% fi = '/scratch/kanga/runs/20180208.GMPAS-IAF.T62_oRRS18v3.anvil/transects/oRRS18to6v3.171116_northatlantic_03.nc';
fi = '/scratch/kanga/runs/20180208.GMPAS-IAF.T62_oRRS18v3.anvil/transects/oRRS18to6v3.171116.Transects_masks.nc';
fi = '/scratch/kanga/runs/EC60to30_G_case/transects/oEC60to30v3_sections_northatlantic_05.nc';
fi = '/scratch/kanga/runs/20180208.GMPAS-IAF.T62_oRRS18v3.anvil/transects/oRRS18to6v3.171116_northatlantic_05.nc';
fi = '/scratch/kanga/runs/oRRS18to6v3.171116_northatlantic_05.nc';

ncread(fi, 'transectNames')';
%%
k = 7;

edge_inds = ncread(fi, 'transectEdgeGlobalIDs', [1,k], [Inf,1]);
edge_inds(edge_inds==0) = [];

i = 1;
%%
verticesOnEdge = ncread(run(i).mesh_fi, 'verticesOnEdge');
vert_inds = verticesOnEdge(:, edge_inds);

lonVertex = rad2deg(ncread(run(i).mesh_fi, 'lonVertex'));
lonVertex(lonVertex>180) = lonVertex(lonVertex>180) - 360;
latVertex = rad2deg(ncread(run(i).mesh_fi, 'latVertex'));
%%
lon = lonVertex(vert_inds(2,:));
lat = latVertex(vert_inds(2,:));


figure
line(lon, lat)


%%
transect(1).name = 'AtlanticCUSP8contour15km';
% transect(1).name = 'AtlanticCrossingat30N';
% transect(2).name = 'AtlanticCrossingat45N';

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

%%
j = 1;

for i = 1:length(run)
    figure
    
    X = repmat(cumsum(run(i).T(j).dx)', [length(run(i).T(j).dz),1]);
    Z = repmat(cumsum(run(i).T(j).dz), [1,length(run(i).T(j).dx)]);
    
    pcolor(X, Z, mean(run(i).T(j).velocity,3) .* run(i).T(j).mask); shading flat

    set(gca,'ydir','reverse')

    crange = [-1 1] * 0.03;
    caxis(crange)
    colormap(cbrewer('div','RdBu',40,'pchip'))
    
end

%%


