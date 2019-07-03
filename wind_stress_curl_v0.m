% Kevin Rosa
% June 6, 2019
%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'CUSP8';
run(i).dir = sprintf('/scratch/kanga/runs/standalone_hoch_mesh/%s/',run(i).name);
run(i).data_fi = fullfile(run(i).dir, 'forcing_data.nc');
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
i = i+1;
run(i).name = 'NA8';
run(i).dir = sprintf('/scratch/kanga/runs/standalone_hoch_mesh/%s/',run(i).name);
run(i).data_fi = fullfile(run(i).dir, 'forcing_data.nc');
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
i = i+1;
run(i).name = 'EC60to30';
run(i).dir = sprintf('/scratch/kanga/runs/standalone_hoch_mesh/%s/',run(i).name);
run(i).data_fi = fullfile(run(i).dir, 'forcing_data.nc');
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');

%%
dx = 0.1;
xrange = [-80 -15];
yrange = [23 37];

lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

for i = 1:length(run)
    
    [run(i).LON, run(i).LAT] = make_lonlat_matrix(lon_vec, lat_vec);

    run(i).mask = compute_mask(run(i).mesh_fi, run(i).LON, run(i).LAT);
    
    [~,~,run(i).u] = mpas_to_lonlat_meshgrid('windStressZonal', run(i).mesh_fi, run(i).data_fi, lon_vec, lat_vec, 1);
    [~,~,run(i).v] = mpas_to_lonlat_meshgrid('windStressMeridional', run(i).mesh_fi, run(i).data_fi, lon_vec, lat_vec, 1);
    
    run(i).u = run(i).u .* run(i).mask;
    run(i).v = run(i).v .* run(i).mask;
end

%%
for i = 1:length(run)
    
    nans = NaN(size(run(i).LON));
    
    u = run(i).u;
    v = run(i).v;
    
    dudx = nans;    dvdy = nans;
    dudx(1:end-1,:) = u(2:end,:) - u(1:end-1,:);
    dvdy(:,1:end-1) = v(:,2:end) - v(:,1:end-1);
    
    run(i).div = dudx + dvdy;
end
    

%%
figure
run_inds = 1:3;
m = 1;
n = length(run_inds);
p = 1;
for i = run_inds
    subplot(m,n,p)
    pcolor(run(i).LON, run(i).LAT, run(i).div); shading flat
    colorbar
    title(run(i).name)
    caxis(1e-3*[-2 4])
    p = p+1;
end

%%
for i = 1:length(run)
    
    nans = NaN(size(run(i).LON));
    
    u = run(i).u;
    v = run(i).v;
    
    dvdx = nans;    dudy = nans;
    dvdx(1:end-1,:) = v(2:end,:) - v(1:end-1,:);
    dudy(:,1:end-1) = u(:,2:end) - u(:,1:end-1);
    
    run(i).curl = dvdx - dudy;
end

%%
figure
run_inds = 1:3;
m = 1;
n = length(run_inds);
p = 1;
for i = run_inds
    subplot(m,n,p)
    pcolor(run(i).LON, run(i).LAT, run(i).curl); shading flat
    colorbar
    title(sprintf('%s: %.3f\n', run(i).name, nansum(run(i).curl(:))))
    caxis(1e-3*[-10 4])
    p = p+1;
end

%%
for i = 1:length(run)
    fprintf('%s: %.3f\n', run(i).name, nansum(run(i).curl(:)))
end
