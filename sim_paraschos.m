clear all; clc;

%% User input
dt = 0.001;
controlGainsFilename = './Data/controlGains_X-tremeFinal.mat';
iiwaUrdfPath = './iiwa_description/urdf/iiwa14.urdf';

%% Init sim
[lbr, q_0, t_end, feedback_gain_K_t, feed_forward_gain_k_t, epsilon_u]...
    = init_sim(controlGainsFilename, iiwaUrdfPath, dt);

%% Sim
disp(strcat('starting at ',{' '}, string(datetime('now'))))
t_start = tic;

simOut = sim('iiwa_paraschos','SimulationMode','normal','StopTime',...
    string(t_end));

toc(t_start)
t_final = toc(t_start);
if t_final > 120
    disp('computer so slooooooooow...')
else
    disp('you so fast, really fast.')
end
disp(strcat('ending at ',{' '}, string(datetime('now'))))

%% Plot
simOut.plot

%% Done
disp('done')


