function [filtereddata,interpolatedtime] = lowpassfilter(time_in,data_in,cutoffperiod_hours)
% LOWPASSFILTER(time_in,data_in,cutoffperiod_hours)
% -input time vector must be in datenum format.
% -default cut-off period is 33 hours
% -interpolates to sampling frequency; NaNs and gaps will be NaNs.
if nargin < 2
    error('don''t forget the time vector')
end

data_in = double(data_in);  % filtfilt input arguments must be 'double'

n = 3; 
dt_days = min(diff(time_in));
dt_s    = dt_days*24*60*60;    % dt (s)
f_s = 1/dt_s;
tidefreq = 1/(33*60*60);
if nargin > 2
    tidefreq = 1/(cutoffperiod_hours*60*60);
end
Wn  = tidefreq/(0.5*f_s); 
[b, a]  = butter(n, Wn, 'low');

T       = time_in(1):dt_days:time_in(end);
data    = interp1(time_in, squeeze(data_in), T);
good    = find(~isnan(data));


filtereddata    = NaN(length(data),1);
filtereddata(good) = filtfilt(b,a, data(good));

interpolatedtime = T;
end

