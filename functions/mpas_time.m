function time = mpas_time(fi, t_ind, time_var_name)
%MPAS_TIME 
% time = mpas_time(fi, t_ind, time_var_name)
%
% 'time_var_name' is optional. Default is 'xtime'.
%   
%
% Kevin Rosa
% June 4, 2019

if ~exist('time_var_name', 'var')
    time_var_name = 'xtime';
end


if strcmp(t_ind,'all')
    time_str = ncread(fi,time_var_name,[1,1],[Inf,Inf]);
    time = datenum(time_str','yyyy-mm-dd_HH:MM:SS');
    
else
    time = NaN(size(t_ind));
    for t = t_ind(:)'
        time_str = ncread(fi,time_var_name,[1,t],[Inf,1]);
        time(t) = datenum(time_str','yyyy-mm-dd_HH:MM:SS');
    end
end
