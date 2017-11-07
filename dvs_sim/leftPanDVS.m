function [ xs, ys, ts, ps ] = leftPanDVS( rate, Tms, size )
% LEFTPANDVS create a series of bars simlating a robot panning left (scene
% moving right). 
%
%   Parameters:
%       rate - scale param for number of events (in units spikes/ms)
%       Tms - Length of simulation (ms)
%
%   Examples:
%   [xs, ys, ts, ps] = leftPanDVS( 2*1e-3, 5000, 128 )

%size = 10;
num_spikes = rate * Tms;
bar_pos = 1;    % Position around which to generate spikes
bar_spacing = ceil(size/4);  % 40% away from eachother
bar_speed = 1000; % in pixels/second
bar_speed_ms = bar_speed / 1000;
noise_factor = 1;%size/8;

ts = sort(ceil(rand(num_spikes, 1) * Tms)*1000);
noise = rand(num_spikes, 1) * noise_factor;
xs = mod(ceil((bar_speed_ms .* (ts / 1000)) + noise), size);
ys = ceil(rand(num_spikes, 1) * size) - 1;
ps = ceil(rand(num_spikes, 1) * 2 - 1);


%plot3(xs, ys, ts, '.')
%xlabel('xs');
%ylabel('ys');

%xs = ceil(rand(num_spikes, 1) * size);
%ys = ceil(rand(num_spikes, 1) * size);
%ts = sort(ceil(rand(num_spikes, 1) * Tms)*1000);
%ps = ceil(rand(num_spikes, 1) * 2 - 1);

end

