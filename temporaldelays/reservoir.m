function [ delays, post ] = reservoir( Ne, Ni, D, M )
%% RESERVOIR - Initialise a random reservoir with delays
%   This is just Izhikevichs code refactored for my new software framework
%   Example useage:
%       [ delays, post ] = reservoir( 800, 200, 20, 100 )

N = Ne + Ni;
delays = cell(N,D);
post = zeros(N,M);
for i=1:Ne
    p=randperm(N);
    post(i,:)=p(1:M);
    for j=1:M
        delays{i, ceil(D*rand)}(end+1) = j;  % Assign random exc delays
    end;
end;
for i=Ne+1:N
    p=randperm(Ne);
    post(i,:)=p(1:M);
    delays{i,1}=1:M;                    % all inh delays are 1 ms.
end;

end