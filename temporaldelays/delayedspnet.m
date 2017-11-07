%% DELAY2LAYER
clear
addpath('../dvs_sim/')
rand('seed',1);

sim_time_ms = 0.038 * 1000;
resolution = 32;
input_size = resolution * resolution;

field_size = 8;
delay_distribution = repmat(1:field_size, field_size, 1);
N_fields_X = resolution / field_size;
N_fields_Y = resolution / field_size;

M = field_size * field_size;
D = 20;
N_inp = input_size;
N_hid = N_fields_X * N_fields_Y; % 1 neuron per field
N = N_inp + N_hid;
a = [0.02*ones(N_inp,1);    0.1*ones(N_hid,1)];
d = [   8*ones(N_inp,1);    2*ones(N_hid,1)];
sm=10;

%% Loading data and preprocessing
rate = 25e0;
[xs, ys, ts, ps] = leftPanDVS( rate, sim_time_ms, resolution );
xs = xs + 1; ys = ys + 1;           % Correct zero indexing
ts = ceil(ts / 1000);               % Convert ts to ms
inp = sub2ind([resolution,resolution], xs, ys);   % Flattern 2D array to neuron indexes
% Now ind(i) spiked at ts(i)
% size(ind) == size(ts)
subplot(2, 2, 1);
plot3(xs, ys, ts, '.k')
title('Input data')
xlabel('xs');
ylabel('ys');
zlabel('ts');


%% Build Model
%[delay_xs, delay_ys] = ind2sub([N_fields_X, N_fields_Y], 1:N_hid);
[delay_xs, delay_ys] = ind2sub([field_size, field_size], 1:M);

delays = cell(N,D);
post = zeros(N, M); 
for j=1:N_hid
    post_idx = j + N_inp;

    % i's jth postsynaptic neuron is post(i, j)
    %post(i,) = ;  % must be in last 200
    %jx = repmat([1:field_size]', field_size, 1);
    %jy = repmat([1:field_size]', field_size, 1);
    %delay_xs, delay_ys = ind2sub([N_fields_X, N_fields_Y], 1:N_hid);
    [jx_off, jy_off] = ind2sub([N_fields_X, N_fields_Y], j);
    jx = (jx_off-1)*field_size + delay_xs;
    jy = (jy_off-1)*field_size + delay_ys;
    
    %[jx; jy]
    inp_idxs = sub2ind([resolution,resolution], jx, jy);
    post(inp_idxs,:) = post_idx; 
    
    for i=1:field_size
        % Delay from i to j (j is index into post), with value D*rand
        %delays{i, ceil(D*rand)}(end+1) = j; 
        % post_j is the indices in post that i connect to.
        post_j = (i-1) * field_size + 1 : i * field_size;
        %delays(inp_idxs(post_j), field_size - i + 1) = num2cell(post_j); 
        delays(inp_idxs(post_j), i ) = num2cell(post_j);
    end;
end

%% Execution
w = 6;
s=[w*ones(N_inp, M);w*ones(N_hid, M)];         % synaptic weights
sd=zeros(N,M);                          % their derivatives

% Make links at postsynaptic targets to the presynaptic weights
pre = cell(N,1);
aux = cell(N,1);
for i=N_inp:N_hid
    for j=1:D
        for k=1:length(delays{i,j})
            pre{post(i, delays{i, j}(k))}(end+1) = N*(delays{i, j}(k)-1)+i;
            aux{post(i, delays{i, j}(k))}(end+1) = N*(D-1-j)+i; % takes into account delay
        end;
    end;
end;
  

STDP = zeros(N,1001+D);
v = -65*ones(N,1);                      % initial values
u = 0.2.*v;                             % initial values
firings=[-D 0];                         % spike timings

v_trace = -65*ones(N, 1000);

disp('Start simulation');
sim_time_s = ceil(sim_time_ms / 1000);
for sec=1:sim_time_s                      
    sec
  for t=1:1000                          % simulation of 1 sec
    I=zeros(N,1); 
    %I(ceil(N*rand))=20;                 % random thalamic input 
    
    % collect input from sensor
    in_fired = zeros(input_size, 1);      % input neurons firing
    time = (sec-1) * 1000 + t; 
    %in_fired(inp(ts <= time & ts > time - dt)) = 1;
    in_fired(inp(ts == time)) = 1;
    %size(inp(ts == time))
    v_trace(:, t) = v;
    
    
    fired = [find(v>=30); inp(ts == time)];                % indices of fired neurons
    %fired = in_fired;  % Input should only fire when added here
    %size(fired)
    v(fired)=-65;  
    u(fired)=u(fired)+d(fired);
    % STDP 
    STDP(fired,t+D)=0.1;
    for k=1:length(fired)
      sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
    end;
    
    % firings is in the form [time, idx; time, idx; ...]
    firings=[firings;t*ones(length(fired),1),fired];
    k=size(firings,1);
    % Starting at most recent firing, go back to when the largest delay
    % could be from
    while firings(k,1)>t-D
      idx = firings(k,2);
      time_fired = firings(k,1);
      delay_post_idxs=delays{idx ,t-time_fired+1};
      post_idxs = post(idx, delay_post_idxs);
      I(post_idxs)=I(post_idxs)+s(idx, delay_post_idxs)';
      sd(idx, delay_post_idxs)=sd(idx,delay_post_idxs)-1.2*STDP(post_idxs,t+D)';
      k=k-1;
    end;
    v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical 
    v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time 
    u=u+a.*(0.2*v-u);                   % step is 0.5 ms
    STDP(:,t+D+1)=0.95*STDP(:,t+D);     % tau = 20 ms
  end;
  %% After 1 second plot how the second layer reacted during time
  subplot(2, 2, 3);
  filter = find(firings(:, 2) > N_inp);
  l2_spike_idxs = firings(filter, 2);
  l2_spike_times = firings(filter, 1);
  [l2_xs, l2_ys] = ind2sub([N_fields_X, N_fields_Y], l2_spike_idxs - N_inp);
  plot3(l2_xs, l2_ys, l2_spike_times, '.k');
  title('Left detector layer')
  xlabel('xs');
  ylabel('ys');
  zlabel('ts');
  
  subplot(2, 2, 2);
  plot(firings(:,1),firings(:,2),'.');
  axis([0 1000 0 N]); drawnow;
  STDP(:,1:D+1)=STDP(:,1001:1001+D);
  post_idxs = find(firings(:,1) > 1001-D);
  firings=[-D 0;firings(post_idxs,1)-1000,firings(post_idxs,2)];
  s(1:N_inp,:)=max(0,min(sm,0.01+s(1:N_inp,:)+sd(1:N_inp,:)));
  sd=0.9*sd;
  

  
end;


