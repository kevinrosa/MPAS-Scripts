% Kevin Rosa
% August 5, 2019

% add AVHRR
% add AVISO

%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = 'EC60to30_G_case';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).color = rgb('blue');
run(i).dx = 60;
i = i+1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).color = rgb('black');
run(i).dx = 8;
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).color = rgb('red');
run(i).dx = 8;
i = i+1;


%% Settings
region = 'cc0';

% download range:
xrange = [-135 -116];
yrange = [30 43];

% which .mat files to read:
target_string = 'timeall_lon-143.0to-116.0_lat_30.0to43.0_highFrequency';

% which fields to read from .mat files:
% FIELDS3D = {'ssh','temperatureAtSurface','barotropicSpeed','relativeVorticityAt250m'};
FIELDS3D = {'ssh','temperatureAtSurface'};

%% Load
for i = 1:length(run)
    D = dir(sprintf('*%s*%s.mat', target_string, run(i).code));
    run(i).fi = D.name;


    % determine spatial indices
    tmp = load(run(i).fi, 'LON','LAT');
    xi  = tmp.LON(:,1)>=xrange(1) & tmp.LON(:,1) <=xrange(2);
    eta = tmp.LAT(1,:)>=yrange(1) & tmp.LAT(1,:) <=yrange(2);

    run(i).LON = tmp.LON(xi,eta);
    run(i).LAT = tmp.LAT(xi,eta);

    % 1D fields
    for F = {'time'}
        tmp = load(run(i).fi, F{1});
        run(i).(F{1}) = tmp.(F{1});
    end

    % additional 2D fields
    for F = {'mask'}
        tmp = load(run(i).fi, F{1});
        run(i).(F{1}) = tmp.(F{1})(xi,eta);
    end

    % 3D fields
    for F = FIELDS3D
        tmp = load(run(i).fi, F{1});
        run(i).(F{1}) = tmp.(F{1})(:,xi,eta);
    end
    
    % pre de-trended ssh
    run(i).ssh_raw = run(i).ssh;
end

%% calculate grid cell widths
var_res_ind = 2;
transition_contour = 12;

[~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;


% %% read AVHRR
% i = 4;  % index in 'run' struct
% 
% data_dir = '/scratch/kanga/AVHRR/';
% 
% D = dir(fullfile(data_dir, 'avhrr-only-v2.year*daily.nc'));
% dt = 5;
% 
% for year = 1:length(D)
%     fi = fullfile(D(year).folder, D(year).name);
%     
%     if year == 1
%         all.lon = ncread(fi, 'lon');
%         all.lat = ncread(fi, 'lat');
%         % convert to -180 180 coords
%         all.lon(all.lon>180) = all.lon(all.lon>180) - 360;
%         
%         xi  = find(all.lon>=xrange(1) & all.lon<=xrange(2));
%         eta = find(all.lat>=yrange(1) & all.lat<=yrange(2));
%         
%         LON = repmat(all.lon(xi), [1, length(eta)]);
%         LAT = repmat(all.lat(eta)', [length(xi), 1]);
%         
%         sst = [];
%     end
%     
%     new = squeeze(ncread(fi, 'sst', [xi(1),eta(1),1,1,1], [length(xi),length(eta),1,1,Inf], [1,1,1,1,dt]));
%     
%     sst = cat(1, permute(new, [3,1,2]), sst);
%     
%     fprintf('%s\n', D(year).name)
% end
%     
% run(i).temperatureAtSurface = sst;
% run(i).LON = double(LON);
% run(i).LAT = double(LAT);
% %%
% run(i).dx = 25;
% run(i).color = rgb('green');

%% read sst
i = 4;
fi = '/scratch/kanga/AVISO/dataset-armor-3d-rep-weekly_1565113697097.nc';
lon = ncread(fi, 'longitude');
lon(lon>180) = lon(lon>180) - 360;
lat = ncread(fi, 'latitude');

run(i).LON = repmat(lon(:), [1,length(lat)]);
run(i).LAT = repmat(lat(:)', [length(lon),1]);

run(i).time = ncread(fi,'time') + datenum(1950,1,1);

run(i).temperatureAtSurface = permute(squeeze(ncread(fi,'to')),[3,1,2]);

run(i).dx = 25;
run(i).color = rgb('green');

%% read sst
i = 5;
fi = '/scratch/kanga/AVISO/dataset-armor-3d-rep-weekly_1565113697097.nc';
lon = ncread(fi, 'longitude');
lon(lon>180) = lon(lon>180) - 360;
lat = ncread(fi, 'latitude');

run(i).LON = repmat(lon(:), [1,length(lat)]);
run(i).LAT = repmat(lat(:)', [length(lon),1]);

run(i).time = ncread(fi,'time') + datenum(1950,1,1);

run(i).ssh = permute(squeeze(ncread(fi,'zo')),[3,1,2]);

run(i).dx = 25;
run(i).color = rgb('green');

%% read AVISO
% i = 5;
% fi = '/scratch/kanga/AVISO/dataset-duacs-rep-global-merged-allsat-phy-l4_1565036020610.nc';
% lon = ncread(fi, 'longitude');
% lon(lon>180) = lon(lon>180) - 360;
% lat = ncread(fi, 'latitude');
% 
% run(i).LON = repmat(lon(:), [1,length(lat)]);
% run(i).LAT = repmat(lat(:)', [length(lon),1]);
% 
% run(i).time = ncread(fi,'time') + datenum(1950,1,1);
% 
% run(i).ssh = permute(ncread(fi,'sla'),[3,1,2]);
% 
% run(i).dx = 25;
% run(i).color = run(i-1).color;

%%


%% Determining track line
% %%
% figure
% m_proj('lambert','long',xrange,'lat',yrange);
% % m_proj('lambert','long',[-140 -70],'lat',[10 60]);
% hold on
% m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
% m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')
% 
% i = var_res_ind;
% m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)
% 
% %%
% set(L,'visible','off')

r = 1300;  % km
center = [-112 41.5];

th = 0:pi/500:2*pi;
xx = r * cos(th);
yy = r * sin(th);

A = 6371 * pi/180;  % dist (km) per degree latitude

lat = center(2) + yy ./ A;
lon = center(1) + xx ./ (A * cosd( (lat+center(2))/2 ));

% discard unwanted indices
lat_range = [31 42];
inds = lat >= lat_range(1) & lat <= lat_range(2) & lon < center(1);

lon = fliplr( lon(inds) );
lat = fliplr( lat(inds) );

% L = m_line(lon,lat,'color','r');


%% calculate trackline in x-y space -- define track line by constant spacing
spacing = 8;  % spacing in km
x = NaN(size(lon));
y = NaN(size(lat));
along_track_dist = NaN(size(lat));
x(1) = 0;
y(1) = 0;
along_track_dist(1) = 0;
for k = 2:length(lon)
    [dx,dy] = lonlat_to_dxdy(lon(k-1), lat(k-1), lon(k), lat(k));
    x(k) = x(k-1) + dx;
    y(k) = y(k-1) + dy;
    along_track_dist(k) = along_track_dist(k-1) + sqrt(dx.^2 + dy.^2);
end

% trackline structure
constant_spacing = 0:spacing:along_track_dist(end);
T.x = interp1(along_track_dist, x, constant_spacing);
T.y = interp1(along_track_dist, y, constant_spacing);

T.lon = interp1(along_track_dist, lon, constant_spacing);
T.lat = interp1(along_track_dist, lat, constant_spacing);


%% check the spacing
for k = 2:length(T.lon)
    [dx,dy] = lonlat_to_dxdy(T.lon(k-1),T.lat(k-1),T.lon(k),T.lat(k));
    fprintf('%.3f km\n', sqrt(dx^2 + dy^2))
end

% sweet!


%% figure showing trackline and transition contour
figure
set(gcf,'color','w','position',[313 496 633 456])

m_proj('lambert','long',[-131 -120],'lat',[31 43]);
hold on
m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid_transparent('box','fancy','tickdir','out','fontsize',12,'linestyle','none')

% i = var_res_ind;
% m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)

m_line(T.lon, T.lat, 'color','r','linewidth',2)

%%
set(gcf,'color','w')
set(gca,'color','none');
export_fig('figures/spectra/trackline_05.png', '-m3','-painters','-transparent','-nocrop')
%
%
%% interpolate data to trackline
for F = {'ssh','temperatureAtSurface'}
field = F{1};
for i = 1:length(run)
    try 
        nt = length(run(i).(field)(:,1,1));
    
        run(i).T.(field) = NaN(nt, length(T.lon));
        for t = 1:nt
            run(i).T.(field)(t,:) = interp2(run(i).LON', run(i).LAT', squeeze(run(i).(field)(t,:,:))', T.lon, T.lat);
        end
    
    catch
        run(i).T.(field) = [];
    end
end
end

%%
run(4).T.ssh = NaN;
run(5).T.temperatureAtSurface = NaN;

%% create one long "timeseries" separated by zeros
% https://www.mathworks.com/matlabcentral/answers/43548-how-to-create-power-spectral-density-from-fft-fourier-transform

pad_dist = 1500; %500; % km

for F = {'temperatureAtSurface', 'ssh'}
    field = F{1};

    figure
    set(gcf,'color','w')
    for i = 1:length(run)

        pad_length = round(pad_dist/spacing);  % put this many zeros between each track's data 
        
        if ~isempty(run(i).T.(field))
        nt = 1:length(run(i).T.(field)(:,1));  % number of tracks
        lt = length(T.lon);  % length of a track

        zero_vec = zeros(pad_length,1);

        data = zero_vec;
        for t = 1:nt
    %         data = cat(1, data, run(i).T.(field)(t,:)' - mean(run(i).T.(field)(t,:)));
    %         data = cat(1, data, detrend(run(i).T.(field)(t,:)'));
            new_vec = run(i).T.(field)(t,:)';
            data = cat(1, data, detrend(new_vec) .* hanning(length(new_vec)));
            data = cat(1, data, zero_vec);
        end
        end

        dx = spacing * 1000;  % m
        Fs = 1/dx; % Sampling frequency 
        L = length(data); % Length of signal 
        tt = (0:L-1)*dx; % Time vector


        xdft = fft(data);
        Pxx = 1/(L*Fs)*abs(xdft(1:floor(length(data)/2+1))).^2;

        wavenumber = 0:Fs/L:Fs/2;

        % remove wavelengths smaller than 2 * dx
        inds = 1./wavenumber < (run(i).dx*2*1000);
        Pxx(inds) = NaN;

        line(log10(wavenumber),log10(Pxx),'color',run(i).color,'linewidth',2);
    end


    xrange = [10 1200] * 1e3;  % wavelength (m)
    xlim( sort(log10(1./xrange)) );


    wavelength = [2500, 1000, 500, 250, 100, 50, 25, 10] * 1e3;

    xticks = log10(1./wavelength);
    set(gca,'xtick',xticks);

    set(gca,'xticklabel',num2str(round(wavelength'* 1e-3)))

    xlabel('Wavelength (km)')

    yticks = get(gca,'ytick');
    yticklabel = {};
    for k = 1:length(yticks)
        yticklabel{k} = sprintf('10^{%i}', yticks(k));
    end
    set(gca,'yticklabel',yticklabel, 'fontsize',12)


    grid on
    
    % SAVE
    set(gca,'color','none');
    export_fig(sprintf('figures/spectra/spectra_hanning_v05_%s.png',field), '-m3','-painters','-transparent')
end


