% Kevin Rosa
% June 14, 2019
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
% i = i+1;
% run(i).name = 'Low-resolution G case';
% run(i).short_name = 'low-res';
% run(i).code = '20180305.GM600.T62_oECv3.eos';
% run(i).mesh_fi = sprintf('/scratch/kanga/runs/%s/mpaso.rst.0050-01-01_00000.nc',run(i).code);
% run(i).color = rgb('blue');

%%
for i = 1:length(run)
    run(i).fi = sprintf('data/transects/%s/transport.floridaCuba.nc',run(i).code);
    
    run(i).time = datenum(ncread(run(i).fi, 'Time'), 0,0); 
    run(i).transport = ncread(run(i).fi, 'Transport');
end

%% mean and std
i = 1;
fprintf('%.2f\n', mean(run(i).transport))
fprintf('%.2f\n', std(run(i).transport))

%%
figure
set(gcf,'color','w')

for i = 1:length(run)
    
    line(run(i).time-run(i).time(1), run(i).transport, 'color',run(i).color,'linewidth',2)
end
xlim([datenum(0,6,1), datenum(22,1,1)])
datetick('x','keeplimits')
set(gca,'fontsize',14)

title('G case Transport through Florida-Cuba')
ylabel('Transport (Sv)')
xlabel('Time (years)')

legend(run(:).short_name)

%%
saveas(gcf,'figures/transport/transport_floridacuba_2runs_v1.png')
