% Kevin Rosa
% July 22, 2019

%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = 'EC60to30_G_case';
run(i).dir = sprintf('/scratch/kanga/runs/%s/',run(i).code);
run(i).data_fi = fullfile(run(i).dir, 'eke_oEC60to30v3_to_0.1x0.1degree/mpaso_ANN_004001_006012_climo.nc');
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).color = rgb('blue');
i = i+1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/',run(i).code);
run(i).data_fi = fullfile(run(i).dir, 'eke_oNAEC60to30cr8L60v1_to_0.1x0.1degree/mpaso_ANN_002301_003712_climo.nc');
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).color = rgb('black');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/',run(i).code);
run(i).data_fi = fullfile(run(i).dir, 'eke_oRRS18to6v3.171116_to_0.1x0.1degree/mpaso_ANN_000201_001912_climo.nc');
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).color = rgb('red');


%% Settings
region = 'cal';

% download range:
xrange = [-131 -120];
yrange = [31 43];

%% read netcdf
FIELDS = {'eke'};%,'timeMonthly_avg_velocityZonal','timeMonthly_avg_velocityMeridional'};
for i = 1:length(run)
    all.lon = ncread(run(i).data_fi, 'lon');
    all.lat = ncread(run(i).data_fi, 'lat');    
    xi = find(all.lon >= xrange(1)-0.05 & all.lon <= xrange(2)+0.05);
    eta = find(all.lat >= yrange(1)-0.05 & all.lat <= yrange(2)+0.05);
    
    for F = FIELDS
        run(i).(F{1}) = ncread(run(i).data_fi, F{1}, [xi(1),eta(1)], [length(xi),length(eta)]);
    end
    run(i).LON = repmat(all.lon(xi),  [1,length(eta)]);
    run(i).LAT = repmat(all.lat(eta)',[length(xi),1]);
end

%% calculate grid cell widths
var_res_ind = 2; 
transition_contour = 12;

i = var_res_ind;
[~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;


%% plotting
w_mmap = 1;
run_inds = 1:3;
FIELD = 'eke';

version_code = 'v1';

crange = [0 220];
dc = 20;
conts = crange(1):dc:crange(2);
% cmap = flipud(cbrewer('div','Spectral',length(conts)-1,'pchip'));
% cmap = flipud(cmocean('thermal',length(conts)-1));
cmap = cmocean('matter',length(conts)-1);
bins = crange(1):dc:crange(2)*2;


for i = run_inds
figure
set(gcf,'color','w','position',[313 496 633 456])

if w_mmap == 0
    contourf(run(i).LON, run(i).LAT, run(i).(FIELD), bins,'linecolor','none')
    set(gca,'xlim',xrange,'ylim',yrange)
    colorbar

else 
    m_proj('lambert','long',xrange,'lat',yrange);
    hold on

    m_contourf(run(i).LON, run(i).LAT, run(i).(FIELD), bins,'linecolor','none')

    % add transition region contour
    if i == var_res_ind
        m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.6*[1,1,1],'linewidth',3)
    end

    m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
    m_grid('box','fancy','tickdir','out','fontsize',12)
    
end

caxis(crange)
colormap(cmap)

save_name = sprintf('figures/mapviews/map_%s_%s_%s_%s.png', region, FIELD, run(i).short_name, version_code);
set(gca,'color','none');
export_fig(save_name, '-m3','-painters','-transparent')
end

%% save one with colorbar
cb = colorbar;
ylabel(cb, 'Surface eddy kinetic energy (cm2 s-2)', 'fontsize',12)

save_name = sprintf('figures/mapviews/map_%s_%s_colorbar_%s.png', region, FIELD, version_code);
set(gca,'color','none');
export_fig(save_name, '-m3','-painters','-transparent')


