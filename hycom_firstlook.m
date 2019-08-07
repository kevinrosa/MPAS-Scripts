fi = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/2000';

t = 1000;
d = 1;

xrange = [-79 -55];
yrange = [35 48];

lon = ncread(fi,'lon');
lat = ncread(fi,'lat');

xi = find(lon>xrange(1) & lon<xrange(2));
eta = find(lat>yrange(1) & lat<yrange(2));

LON = repmat(lon(xi),[1,length(eta)]);
LAT = repmat(lat(eta)',[length(xi),1]);

sst = ncread(fi,'water_temp',[xi(1),eta(1),d,t],[length(xi),length(eta),1,1]);
% ssh = ncread(fi,'surf_el',[xi(1),eta(1),t],[length(xi),length(eta),1]);

time = ncread(fi,'time',t,1) + datenum(2000,1,1);

%% bathymetry
% ftp://ftp.hycom.org/datasets/GLBb0.08/expt_93.0/topo/
fi = '/scratch/kanga/depth_GLBb0.08_09m11.nc';

lon = ncread(fi,'Longitude');
lat = ncread(fi,'Latitude');
h = ncread(fi, 'depth');

%%
figure
pcolor(LON, LAT, sst); shading flat

cmap = flipud(cbrewer('div','Spectral',60,'pchip'));
colormap(cmap)

title(sprintf('HYCOM reanalysis 3.1 SST %s',datestr(time,'yyyy-mm-dd')))
set(gca,'fontsize',16)

caxis([5 25])
colorbar

%%
hold on
contour(lon, lat, h, 300*[1,1],'color','k','linewidth',2)

grid on

%%
set(gcf,'color','w')
export_fig('figures/hycom_sst_00','-m2')
