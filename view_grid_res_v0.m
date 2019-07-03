addpath(genpath('.'))

i = 1;
run(i).name = 'B case (fully coupled)';
run(i).short_name = 'B-case';
run(i).dir = '/scratch/kanga/runs/A_WCYC1850_ne30_oNAEC60to30cr8L60v1_anvil01/mpaso.hist.am.highFrequencyOutput/';
run(i).mesh_fi = '/scratch/kanga/runs/A_WCYC1850_ne30_oNAEC60to30cr8L60v1_anvil01/mpaso.rst.0001-01-06_00000.nc';
run(i).color = rgb('grey');
i = i+1;
run(i).name = 'G case (CORE atmosphere)';
run(i).short_name = 'G-case';
run(i).dir = '/scratch/kanga/runs/GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01/mpaso.hist.am.highFrequencyOutput/';
run(i).mesh_fi = '/scratch/kanga/runs/GMPAS-IAF_T62_oNAEC60to30cr8L60v1_anvil01/mpaso.rst.0001-08-01_00000.nc';
run(i).color = rgb('blue purple');

%%
for i = 2%:length(run)
    data_fi = run(i).mesh_fi;
    
    dx = 0.1;
    lon_vec = -84:dx:-30;
    lat_vec = 20:dx:60;
        
    t_ind = 0;
    

    [LON, LAT, areaCell] = mpas_to_lonlat_meshgrid('areaCell', run(i).mesh_fi, data_fi, lon_vec, lat_vec, t_ind);
    
end

%%
DX = 2 * sqrt(areaCell ./ pi) * 1e-3;

%%
figure
pcolor(LON, LAT, DX); shading flat
