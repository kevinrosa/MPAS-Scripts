% Kevin Rosa
% July 30, 2019

%%
addpath(genpath('.'))

%%
type = 'transects';
mask_version = 'northatlantic_06';
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
transect(1).name = 'AtlanticCrossingat26.5N';
transect(2).name = 'DavisStrait';
transect(3).name = 'LancasterSound';
transect(4).name = 'NaresStrait';
transect(5).name = 'BarentsSeaOpening';



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
m = load('functions/m_map/private/m_coasts.mat');
%%
figure
line(m.ncst(:,1), m.ncst(:,2),'color','k')
grid on

i = 2;

for j = 1:length(transect)
    line(run(i).T(j).lon, run(i).T(j).lat, 'color','r','linewidth',3)
end

