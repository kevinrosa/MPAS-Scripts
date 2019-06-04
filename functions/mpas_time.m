function time = mpas_time(fi, t_ind)
%MPAS_TIME 
% time = mpas_time(fi, t_ind)
%   
%
% Kevin Rosa
% June 4, 2019

time = NaN(size(t_ind));

for t = t_ind(:)'
    time_str = ncread(fi,'xtime',[1,t],[Inf,1]);
    time(t) = datenum(time_str','yyyy-mm-dd_HH:MM:SS');
end

