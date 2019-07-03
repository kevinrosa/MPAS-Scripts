
% Kevin Rosa
% June 24, 2019
%%
addpath(genpath('.'))

%%
i = 1;
run(i).name = 'Coastally-refined G case';
run(i).short_name = 'var-res';
run(i).code = 'GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01';
run(i).dir = sprintf('/scratch/kanga/runs/%s/',run(i).code);
run(i).fi  = fullfile(run(i).dir, 'mpaso_ANN_002301_003712_climo.nc');
run(i).color = rgb('red');
i = i+1;
run(i).name = 'High-resolution G case';
run(i).short_name = 'high-res';
run(i).code = '20180208.GMPAS-IAF.T62_oRRS18v3.anvil';
run(i).dir = sprintf('/scratch/kanga/runs/%s/',run(i).code);
run(i).fi  = fullfile(run(i).dir, 'mpaso_ANN_000201_001912_climo.nc');
run(i).color = rgb('black');
i = i+1;
run(i).name = 'Low-resolution G case';
run(i).short_name = 'low-res';
run(i).code = '20180305.GM600.T62_oECv3.eos';
run(i).dir = sprintf('/scratch/kanga/runs/%s/',run(i).code);
run(i).color = rgb('blue');

%%
FIELDS = {'lat','lon','eke','timeMonthly_avg_velocityZonal','timeMonthly_avg_velocityMeridional'};
for i = 1:2%length(run)
    for F = FIELDS
        run(i).(F{1}) = ncread(run(i).fi, F{1});
    end
    run(i).LON = repmat(run(i).lon(:), [1,length(run(i).lat)]);
    run(i).LAT = repmat(run(i).lat(:)',[length(run(i).lon),1]);
end

%%
region = 'cc0';
xrange = [-131 -120];
yrange = [31 43];

region = 'gs0';
xrange = [-84 -63];
yrange = [22 38];

region = 'gs2';
xrange = [-77 -67];
yrange = [34 42];



%% plotting
w_mmap = 1;
run_inds = [3,1,2];
FIELD = 'eke';

m = 1;
n = length(run_inds);

version_code = 'v0';

var_res_ind = 1;
transition_contour = 15;
    
figure
set(gcf,'name',FIELD,'position',[13 449 1629 503],'color','w')
pcounter = 1;


if strcmp(region,'cc0')
    crange = [0 250];
    dc = 10;
elseif strcmp(region,'gs0')
    crange = [0 500];
    dc = 20;
elseif strcmp(region,'gs2')
    crange = [0 1500];
    dc = 20;
end


bins = crange(1):dc:crange(2);
cmap = flipud(cbrewer('div','Spectral',length(bins)-1,'pchip'));

for i = run_inds
subplot(m,n,pcounter)

TITLE = sprintf('%s %s', run(i).short_name, FIELD);

if ~isempty(run(i).LON)
    % only bother plotting indices in target lonlat range
    dx = 5;
    xi  = run(i).LON(:,1)>xrange(1)-dx & run(i).LON(:,1)<xrange(2)+dx;
    eta = run(i).LAT(1,:)>yrange(1)-dx & run(i).LAT(1,:)<yrange(2)+dx;

    data = run(i).(FIELD)(xi,eta);
    LON = run(i).LON(xi,eta);
    LAT = run(i).LAT(xi,eta);
    
else
    LON = xrange;
    LAT = yrange;
    data = NaN(2);

end


if w_mmap == 0
    contourf(LON, LAT, data, bins, 'linestyle','none')
    
    set(gca,'xlim',xrange,'ylim',yrange)

else 
    m_proj('lambert','long',xrange,'lat',yrange);
    hold on

    m_contourf(LON, LAT, data, bins,'linecolor','none')

%     % add transition region contour
%     if i == var_res_ind
%         m_contour(run(i).LON, run(i).LAT, run(i).widthCell, transition_contour*[1,1],'linecolor',0.3*[1,1,1],'linewidth',2)
%     end

    m_gshhs_i('patch',0.8*[1,1,1],'edgecolor','k');
    m_grid('linewi',2,'tickdir','out','fontsize',10,'linestyle','none')    
    
end

caxis(crange)
colormap(cmap)
cb = colorbar;



ylabel(cb,'EKE')

title(TITLE)

pcounter = pcounter+1;
end

save_name = sprintf('figures/eke/eke_%s_v0.png', region);
saveas(gcf, save_name)

%%
figure
pcolor(run(i).LON, run(i).LAT, run(i).timeMonthly_avg_velocityMeridional); shading flat

%%
caxis(0.2*[-1 1])
colormap(flipud(cmocean('balance',40)))
