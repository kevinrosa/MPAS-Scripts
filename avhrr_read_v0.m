% Kevin Rosa
% July 11, 2019

%%
addpath(genpath('.'))

%%
data_dir = '/scratch/kanga/AVHRR/';

fi = fullfile(data_dir, 'avhrr-only-v2.year2004daily.nc');

xrange = [-135 -116];
yrange = [30 43];

%%
all.lon = ncread(fi, 'lon');
all.lat = ncread(fi, 'lat');

% convert to -180 180 coords
all.lon(all.lon>180) = all.lon(all.lon>180) - 360;

xi  = find(all.lon>=xrange(1) & all.lon<=xrange(2));
eta = find(all.lat>=yrange(1) & all.lat<=yrange(2));

%%
LON = repmat(all.lon(xi), [1, length(eta)]);
LAT = repmat(all.lat(eta)', [length(xi), 1]);

%%
% sst = ncread(fi, 'sst', [xi(1),eta(1),1], [length(xi),length(eta),Inf]);
sst = squeeze(ncread(fi, 'sst', [xi(1),eta(1),1,1,1], [length(xi),length(eta),1,1,Inf]));

sst = permute(sst, [3,1,2]);

%%
figure
t = 180;
pcolor(LON, LAT, sst(:,:,t+5)); shading flat 

