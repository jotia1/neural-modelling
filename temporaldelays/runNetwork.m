function [ firings ] = runNetwork( delays, post, stimuli )
%% RUNNETWORK - Execute the given model and get firing times of neurons
%   
%   Parameters:
%       delays - delays of each neurons in network
%       post - postsynaptic connections of each neuron in the network
%       stimuli - matrix of [xs, ys, ts, ps]
%
%   Example usage:
%       # TODO - needs lots of params

% TODO 
%   - Need to define an interface for network params such as N_inp
%   - Need to disable STDP for now - later it can be optional


rand('seed',1);

sim_time_ms = 100 * 1000;
%resolution = 32;
%input_size = resolution;
%input_size = 800;

%field_size = 8;
%delay_distribution = repmat(1:field_size, field_size, 1);
delay_distrib = 1:3;
%N_fields_X = resolution / field_size;
%N_fields_Y = resolution / field_size;

%M = 100;
M = size(post, 2);
D = size(delays, 2);
N_inp = 800; % TODO need to get these params some how
N_hid = 200; %Tinput_size - 2; % Lose one on each end
N = N_inp + N_hid;
a = [0.02*ones(N_inp,1);    0.1*ones(N_hid,1)];
d = [   8*ones(N_inp,1);    2*ones(N_hid,1)];
sm=10;
v_thres = -55;

%% Loading data and preprocessing
%[ xs, ys, ts, ps ] = rightDot1D( sim_time_ms, input_size, 1 );
%[ xs, ys, ts, ps ] = twoDot1D( sim_time_ms/2, input_size, 1 );
xs = stimuli(:, 1); ys = stimuli(:, 2); 
ts = stimuli(:, 3); ps = stimuli(:, 4);
xs = xs + 1; ys = ys + 1;           % Correct zero indexing
ts = ceil(ts / 1000);               % Convert ts to ms
inp = xs;  % xs is which neuron to fire

disp('Dont forget to plot xs to make sure it is correctly extracted');




%% Build Model
%[delays, post] = rightDetector(N_inp, N_hid, D);
%[ delays, post ] = reservoir( N_inp, N_hid, D, M );


%% Execution
w = 6;
s=[w*ones(N_inp, M);-5*ones(N_hid, M)];         % synaptic weights
sd=zeros(N,M);                          % their derivatives

% Make links at postsynaptic targets to the presynaptic weights
pre = cell(N,1);
aux = cell(N,1);
for i=N_inp
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
  %  sec
  for t=1:1000                          % simulation of 1 sec
    I=zeros(N,1); 
    
    % collect input from sensor
    %in_fired = zeros(input_size, 1);      % input neurons firing
    time = (sec-1) * 1000 + t; 
    %in_fired(inp(ts <= time & ts > time - dt)) = 1;
    %in_fired(inp(ts == time)) = 0;
    %size(inp(ts == time))
    v_trace(:, t) = v;
    
    
    fired = [find(v>=v_thres); inp(ts == time)'];                % indices of fired neurons
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
  clf
  subplot(2, 1, 1);
  filter = find(firings(:, 2)>N_inp);
  l2_spike_idxs = firings(filter, 2);
  l2_spike_times = firings(filter, 1);
  plot(firings(:,1),firings(:,2),'.','MarkerSize', 8);  % Plot all
  hold on
  plot(l2_spike_times, l2_spike_idxs, '.r','MarkerSize', 8)  % Replot layer 2 in red
  axis([0 1000 0 N]); drawnow;
  title('Full network response')
  xlabel('Time [ms]')
  ylabel('Neuron number')
  legend({'L1', 'L2'})

  
%   subplot(2, 2, 3);
%   plot(firings(:,1),firings(:,2),'.');
%   axis([0 sim_time_ms 0 N]); drawnow;
%   %STDP(:,1:D+1)=STDP(:,1001:1001+D);
%   post_idxs = find(firings(:,1) > 1001-D);
  firings=[-D 0;firings(post_idxs,1)-1000,firings(post_idxs,2)];
  s(1:N_inp,:)=max(0,min(sm,0.01+s(1:N_inp,:)+sd(1:N_inp,:)));
  sd=0.9*sd;
  
  out_neuron = 30;
  in_neurons = [2, 3, 4] + (out_neuron-19);
  subplot(2, 1, 2);
  plot(v_trace(out_neuron, :), 'r')
  hold on
  plot(v_trace(in_neurons(1), :), 'g')
  plot(v_trace(in_neurons(2), :), 'b')
  plot(v_trace(in_neurons(3), :), 'y')
  %refline([0 v_thres]);
  plot( get( gca, 'Xlim' ), [v_thres v_thres], '--k', 'LineWidth',2)
  legend(cellstr(num2str([out_neuron, in_neurons]', '%-d')))
  axis([-Inf 1000 -Inf Inf]); drawnow;
  title('Voltage trace of some dots')
  xlabel('Time [ms]')
  ylabel('Voltage [mV]')
  waitforbuttonpress;
  
end;


end