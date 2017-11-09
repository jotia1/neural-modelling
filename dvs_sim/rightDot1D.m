function [ xs, ys, ts, ps ] = rightDot1D( Tms, ret_length, speed )
%% RIGHTDOT1D - Create a 1D retina with a dot moving with constant velocity
%   
%   Parameters:
%       Tms - Time in ms of the sample to create
%       size - Length of the 1D array (number of pixels in retina)
%       speed - in units of pixels per ms that the dot should move
%
%   Example usage:
%       [ xs, ys, ts, ps ] = rightDot1D( 500, 32, 0.2 )

%dot_speed = 0.2;       % pixels / ms

num_spikes = Tms * speed;
ts = (1:num_spikes) / speed * 1000; % convert to us
xs = mod(0:num_spikes - 1, ret_length);
ys = ones(size(xs)) - 1;  % Make zero indexing to reflect DVS
ps = ceil(rand(size(xs)) * 2 - 1);

end