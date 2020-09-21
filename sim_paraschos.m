clear all; clc;

%% User input
dt = 0.001;
tmp1 = load('controlGains_X-tremeFinal.mat');

%%
t_end = length(tmp1.q_mean)*dt;
q_0 = tmp1.q_mean(1,:);

%%
lbr = importrobot('./iiwa_description/urdf/iiwa14.urdf');
lbr.DataFormat = 'column';
lbr.Gravity = [0, 0, -9.80665];

%% load ProMP

nr_steps = size(tmp1.K,3);
tmp2 = 0:dt:nr_steps*dt-dt;

% feedback_gain_K_t = timeseries(K,nr_steps);
feedback_gain_K_t.time=tmp2;
feedback_gain_K_t.signals.values=tmp1.K;
feedback_gain_K_t.signals.dimensions=size(tmp1.K,[1,2]);

% feed_forward_gain_k_t = timeseries(k,nr_steps);
feed_forward_gain_k_t.time=tmp2;
feed_forward_gain_k_t.signals.values=tmp1.k';
feed_forward_gain_k_t.signals.dimensions=7;

% epsilon_u = timeseries(k,nr_steps);
epsilon_u.time = tmp2;
epsilon_u.signals.values = tmp1.epsilon_u';
epsilon_u.signals.dimensions = 7;

clear tmp1 tmp2
%%

disp(strcat('starting at ',{' '}, string(datetime('now'))))
t_start = tic;
simOut = sim('iiwa_paraschos','SimulationMode','normal','StopTime',string(t_end));
toc(t_start)
t_final = toc(t_start);
if t_final > 120
    disp('computer so slooooooooow...')
end
disp(strcat('ending at ',{' '}, string(datetime('now'))))

disp('done')

%%
simOut.plot

