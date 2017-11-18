
function [  xs, ys, ts, ps ] = spatial(  Tms, ms_per_flash, pattern  )
%SPATIAL Generate spatial data
%   The idea of this function is to generate data that is spatial
%   interesting but temporally boring such that we can exaimine questions
%   of spatial prediction without temporal influence. Data comes back as
%   bars along a 3x3 retina.
%
%   Params
%       Tms - Length of data sent back (ms)
%       ms_per_flash - Frequency of flashes (ms per flash)
%       pattern - the neurons to flash
%
%   Example usage:
%       [  xs, ys, ts, ps ] = spatial(  500, 100, 1:3  )

active = pattern; %active = 1:3;
num_flashs = floor(Tms / ms_per_flash);
xs = repmat(active, 1, num_flashs);
ys = ones(size(xs));
pattern_size = numel(active);
ts = (floor((0:numel(xs)-1)/pattern_size)+1) * ms_per_flash;
ps = ceil(rand(size(xs)) * 2 - 1);

end

