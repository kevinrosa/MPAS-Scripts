% Kevin Rosa
% June 27, 2019

fi = '/scratch/kanga/runs/20180305.GM600.T62_oECv3.eos/mpaso.hist.am.timeSeriesStatsMonthly/mpaso.hist.am.timeSeriesStatsMonthly.0001-08-01.nc';
mesh_fi = '/scratch/kanga/runs/20180305.GM600.T62_oECv3.eos/mpaso.rst.0050-01-01_00000.nc';

%%
%     timeMonthly_avg_velocityMeridional                                                            
%            Size:       60x235160x1
%            Dimensions: nVertLevels,nCells,Time
%            Datatype:   double
%            Attributes:
%                        units     = 'm s^{-1}'
%                        long_name = 'component of horizontal velocity in the northward direction'
%     timeMonthly_avg_velocityZonal                                                                 
%            Size:       60x235160x1
%            Dimensions: nVertLevels,nCells,Time
%            Datatype:   double
%            Attributes:
%                        units     = 'm s^{-1}'
%                        long_name = 'component of horizontal velocity in the eastward direction'
%     timeMonthly_avg_layerThickness                                                                
%            Size:       60x235160x1
%            Dimensions: nVertLevels,nCells,Time
%            Datatype:   double
%            Attributes:
%                        units     = 'm'
%                        long_name = 'layer thickness'
% 

%{
timeMonthly_avg_layerThickness is 0 for levels below bottom.
so can just take column depth to be sum of layerThickness.

add a new variable waterColumnThickness. 
%}

%%
u = ncread(fi,'timeMonthly_avg_velocityMeridional');
v = ncread(fi,'timeMonthly_avg_velocityZonal');
dz = ncread(fi, 'timeMonthly_avg_layerThickness');

%%
ubar = sum(u .* dz) ./ sum(dz);
vbar = sum(v .* dz) ./ sum(dz);



