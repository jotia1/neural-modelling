function [ delays, post ] = rightDetector( N_inp, N_hid, D )

N = N_inp + N_hid;
M = 3;
delays = cell(N, D);
post = zeros(N, M);
delay_pattern = M:-1:1;

%% Do main body of delays
for pre_idx = D : N_inp - D + 1
    delays(pre_idx, 1:D) = num2cell(1:M);  % All middle have 3 synapses
    post(pre_idx, 1:M) = pre_idx + N_inp + [1:3] - M;
end

delays(1, M) = num2cell(1:1); % 1 has one synapse, with delay M
post(1, 1) = 1 + N_inp + [3:3] - M;
delays(2, M-1:M) = num2cell(1:2); % 2 has two synapses
post(2, 1:2) = 2 + N_inp + [2:3] - M;
delays(N_inp, 1) = num2cell(1:1); % N_inp'th has one synapse
post(N_inp, 1) = N_inp + N_inp + [1:1] - M;
delays(N_inp-1, 1:2) = num2cell(1:2); % N_inp-1'th has 2 synpases
post(N_inp-1, 1:2) = N_inp - 1 + N_inp + [1:2] - M;


end