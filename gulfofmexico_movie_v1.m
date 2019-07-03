% Kevin Rosa
% July 2, 2019
%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'CUSP8';
run(i).short_name = run(i).name;
run(i).code = run(i).name;
run(i).dir = sprintf('/scratch/kanga/runs/standalone_hoch_mesh/%s/',run(i).code);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).years = 4:10;
i = i+1;
run(i).name = 'CUSP12';
run(i).short_name = run(i).name;
run(i).code = run(i).name;
run(i).dir = sprintf('/scratch/kanga/runs/standalone_hoch_mesh/%s/',run(i).code);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).years = 4:10;
i = i+1;
run(i).name = 'NA8';
run(i).short_name = run(i).name;
run(i).code = run(i).name;
run(i).dir = sprintf('/scratch/kanga/runs/standalone_hoch_mesh/%s/',run(i).code);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).years = 4:10;
i = i+1;
run(i).name = 'CUSP20';
run(i).short_name = run(i).name;
run(i).code = run(i).name;
run(i).dir = sprintf('/scratch/kanga/runs/standalone_hoch_mesh/%s/',run(i).code);
run(i).mesh_fi = fullfile(run(i).dir, 'init.nc');
run(i).years = 4:10;

%% Settings
% download range
xrange = [-97 -80];
yrange = [18 30];
dx = 0.1;
lon_vec = xrange(1):dx:xrange(2);
lat_vec = yrange(1):dx:yrange(2);

FIELDS = {'ssh'};

for i = 4%1:length(run)
    files = {};
    for year = run(i).years
        dd = dir(fullfile(run(i).dir, sprintf('mpaso.hist.am.highFrequencyOutput.%04i-*',year)));
        files = cat(1, files, fullfile({dd(:).folder}, {dd(:).name})');
    end

    t_ind = 1;
    t_length = length(files);  % number of time indices 

    nans1d = NaN(t_length, 1);
    run(i).time = nans1d;
    
    nans3d = NaN([t_length, length(lon_vec), length(lat_vec)]);
    for F = FIELDS
        run(i).(F{1}) = nans3d;
    end

    tt = 1;
    for m = 1:length(files)
        
        data_fi = files{m};
        
        run(i).time(tt) = mpas_time(data_fi, t_ind);            

        
        % Calculate LON LAT and land-sea mask (only do this once)
        if tt == 1
            [LON, LAT] = make_lonlat_matrix(lon_vec, lat_vec);
            run(i).LON = LON;
            run(i).LAT = LAT;
            
            run(i).mask = compute_mask(run(i).mesh_fi, LON, LAT);
                        
        end
        
        % read data fields
        for F = FIELDS
            [~,~,data] = mpas_to_lonlat_meshgrid(F{1}, run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
            run(i).(F{1})(tt,:,:) = data .* run(i).mask;
        end
        
        % pre-corrected ssh
        run(i).ssh_raw = run(i).ssh;
        
        fprintf('%s %.1f%s\n', run(i).short_name, 100*tt/t_length, '%')
        tt = tt+1;
    end
end


%% de-trend SSH
for i = 1:length(run)
    run(i).ssh_trend = NaN(length(run(i).time),1);
    for t = 1:length(run(i).time)
        data = run(i).ssh_raw(t,:,:);
        run(i).ssh_trend(t) = nanmean(data(:));
    end
    
    run(i).ssh = run(i).ssh_raw - repmat(run(i).ssh_trend,[1,size(run(i).LON)]);
end

%%
crange = [-1 1]*0.4;
dc = 0.01;
bins = 2*crange(1):dc:crange(2)*2;

cmap = flipud(cbrewer('div','Spectral',length(bins)-1,'pchip'));

%%
for i = 4%2%1:length(run)

nFrames = length(run(i).time);
vidObj = VideoWriter(sprintf('mov_loopcurrent_%s_v0.avi',run(i).code));
vidObj.Quality = 100;
vidObj.FrameRate = 12;
open(vidObj)

figure
set(gcf,'color','w')

m_proj('lambert','long',xrange,'lat',yrange);
hold on
m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')

colormap(cmap)

for t = 1:length(run(i).time)
    [~,H] = m_contourf(run(i).LON, run(i).LAT, squeeze(run(i).ssh(t,:,:)), bins,'linecolor','none');
    
    caxis(crange)
    title(datestr(run(i).time(t)))
    
    writeVideo(vidObj, getframe(gcf));
    set(H,'visible','off')
end

close(vidObj);

end
%% 
