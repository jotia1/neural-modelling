function [ xs, ys, ts, ps ] = poissonDVS( rate, Tms )
%FAKEDVS Make a fake DVS signal in a 10x10 square
%   DVS signal is just poisson noise
%
%   Examples:
%   [xs, ys, ts, ps] = poissonDVS( 2*1e-3, 5000 )

size = 10;
num_spikes = rate * Tms * size * size;

xs = ceil(rand(num_spikes, 1) * size);
ys = ceil(rand(num_spikes, 1) * size);
ts = sort(ceil(rand(num_spikes, 1) * Tms)*1000);
ps = ceil(rand(num_spikes, 1) * 2 - 1);

end

