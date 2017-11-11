
function [  xs, ys, ts, ps ] = twoDot1D(  Tms, ret_length, speed  )
%TWODOT1D Summary of this function goes here
%   Detailed explanation goes here

[ xs, ys, ts, ps ] = rightDot1D( Tms, ret_length, speed );

% Flip and duplicate xs
xs_flipped = ret_length - xs - 1;  % flip
xss = [xs; xs_flipped];
xs_two_dot = reshape(xss, 1, numel(xss));
xs = xs_two_dot;

% duplicate others
dup = @(m) reshape(repmat(m, 2, 1), 1, numel(m) * 2);
ys = dup(ys);
ps = dup(ps);
ts = dup(ts);
%[ys, ps, ts] = cell2mat(cellfun(dup, {ys, ps, ts}));
% ys = reshape(repmat(ys, 2, 1), 1, numel(ys) * 2);
% ps = reshape(repmat(ps, 2, 1), 1, numel(ps) * 2);
% ts = reshape(repmat(ts, 2, 1), 1, numel(ts) * 2);

end

