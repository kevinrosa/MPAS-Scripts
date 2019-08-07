% Kevin Rosa
% July 12, 2019
%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0001-08-01_00000.nc',run(i).code);
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/oRRS18to6v3.171116-7.nc',run(i).code);
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = 'EC60to30_G_case';
run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
run(i).color = rgb('blue');

% matObj = matfile('Pacific_years20to36_time1_lon-143.0to-116.0_lat_30.0to43.0_highFrequency_20180208.GMPAS-IAF.T62_oRRS18v3.anvil.mat');
% varlist = who(matObj)

%% Settings
region = 'cc0';

% download range:
xrange = [-135 -116];
yrange = [30 43];

% which .mat files to read:
target_string = 'time1_lon-143.0to-116.0_lat_30.0to43.0_highFrequency';

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

%% de-trend SSH
% for i = 1:length(run)
%     run(i).ssh_trend = NaN(length(run(i).time),1);
%     for t = 1:length(run(i).time)
%         data = run(i).ssh_raw(t,:,:);
%         run(i).ssh_trend(t) = nanmean(data(:));
%     end
%     
%     run(i).ssh = run(i).ssh_raw - repmat(run(i).ssh_trend,[1,size(run(i).LON)]);
% end
% 
% %% pick a time for snapshots
% target_date = datenum(0,07,1);
% for i = 1:length(run)
%     dv = datevec(run(i).time);
%     day_of_year = datenum(0,dv(:,2),dv(:,3));
%     
%     close_times = find( abs(day_of_year-target_date) < 6 );
%     
%     % pick a date in the middle of the run
%     run(i).snap_t_ind = close_times(ceil(length(close_times)/2));
% end



%% calculate grid cell widths
for i = 1:length(run)
    [~, ~, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, run(i).mesh_fi, run(i).LON(:,1), run(i).LAT(1,:), 0);
    run(i).widthCell = 2 * sqrt(areaCell ./ pi) * 1e-3;
end

var_res_ind = 1;
transition_contour = 15;

%% Determining track line
%%
figure
m_proj('lambert','long',xrange,'lat',yrange);
% m_proj('lambert','long',[-140 -70],'lat',[10 60]);
hold on
m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')

i = var_res_ind;
m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)

%%
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

%%
m_line(T.lon,T.lat,'marker','+','color','b')

%% check the spacing
for k = 2:length(T.lon)
    [dx,dy] = lonlat_to_dxdy(T.lon(k-1),T.lat(k-1),T.lon(k),T.lat(k));
    fprintf('%.3f km\n', sqrt(dx^2 + dy^2))
end

% sweet!

%% figure showing trackline and transition contour
figure
m_proj('lambert','long',xrange,'lat',yrange);
hold on
m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')

i = var_res_ind;
m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)

m_line(T.lon, T.lat, 'color','r','linewidth',2)
%%
set(gcf,'color','w')
export_fig('figures/spectra/trackline_00','-m2')

%%
%
%
%% interpolate data to trackline
for F = {'ssh','temperatureAtSurface'}
field = F{1};
for i = 1:length(run)
    run(i).T.(field) = NaN(length(run(i).time), length(T.lon));
    for t = 1:length(run(i).time)
        run(i).T.(field)(t,:) = interp2(run(i).LON', run(i).LAT', squeeze(run(i).(field)(t,:,:))', T.lon, T.lat);
    end
end
end

%% avhrr
sat.temperatureAtSurface = sst;
sat.LON = LON;
sat.LAT = LAT;
field = 'temperatureAtSurface';
sat.T.(field) = NaN(length(sst(:,1,1)), length(T.lon));
for t = 1:length(sst(:,1,1))
    sat.T.(field)(t,:) = interp2(sat.LON', sat.LAT', squeeze(sat.(field)(t,:,:))', T.lon, T.lat);
end


%% create one long "timeseries" separated by zeros
% https://www.mathworks.com/matlabcentral/answers/43548-how-to-create-power-spectral-density-from-fft-fourier-transform

i = 4;
run(i) = sat;
run(i).color = 'k';

pad_dist = 500; % km
field = 'temperatureAtSurface';

figure
for i = 4%1:length(run)
    
    pad_length = round(pad_dist/spacing);  % put this many zeros between each track's data 
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
    
    dx = spacing * 1000;  % m
    Fs = 1/dx; % Sampling frequency 
    L = length(data); % Length of signal 
    tt = (0:L-1)*dx; % Time vector


    xdft = fft(data);
    Pxx = 1/(L*Fs)*abs(xdft(1:floor(length(data)/2+1))).^2;

    wavenumber = 0:Fs/L:Fs/2;
    line(log10(wavenumber),log10(Pxx),'color',run(i).color,'linewidth',2);

end



%% create one long "timeseries" separated by zeros
% https://www.mathworks.com/matlabcentral/answers/43548-how-to-create-power-spectral-density-from-fft-fourier-transform

pad_dist = 500; % km
field = 'temperatureAtSurface';

figure
for i = 1:length(run)
    
    pad_length = round(pad_dist/spacing);  % put this many zeros between each track's data 
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
    
    dx = spacing * 1000;  % m
    Fs = 1/dx; % Sampling frequency 
    L = length(data); % Length of signal 
    tt = (0:L-1)*dx; % Time vector


    xdft = fft(data);
    Pxx = 1/(L*Fs)*abs(xdft(1:floor(length(data)/2+1))).^2;

    wavenumber = 0:Fs/L:Fs/2;
    line(log10(wavenumber),log10(Pxx),'color',run(i).color,'linewidth',2);

end
%%
xrange = [10 1200] * 1e3;  % wavelength (m)
xlim( sort(log10(1./xrange)) );


% xticks = get(gca,'xtick');
% wavelength = 1 ./ (10.^xticks);
% 
% set(gca,'xticklabel',num2str(round(wavelength'./1000)))

% wavelength = [2500, 1000, 250, 100, 25, 10] * 1e3;
% wavelength = [3000:-1000:1000, 500:-100:100, 50:-10:10] * 1e3;
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
set(gca,'yticklabel',yticklabel)


grid on

% saveas(gcf, sprintf('figures/spectra/spectra_hanning_v01_%s.png',field))

%%
figure
t = 1;
field = 'temperatureAtSurface';
for i = 1:length(run)
    line(constant_spacing, detrend(run(i).T.(field)(t,:)), 'color',run(i).color,'linewidth',2);
end

%%

