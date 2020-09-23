function [...
    lbr, q_0, t_end, feedback_gain_K_t, feed_forward_gain_k_t,...
    epsilon_u] ...
    = init_sim(...
    controlGainsFilename, iiwaUrdfPath, dt...
    )
%INIT_SIM Summary of this function goes here
%   Detailed explanation goes here

%% Check args
if nargin < 2
    errID = 'init_sim:BadInput';
    errMsg = 'Number of arguments incorrect. nargin was %d';
    baseException = MException(errID, errMsg, nargin);
    
    causeException = MException(errID, errMsg, nargin);
    baseException = addCause(baseException,causeException);
    
    throw(baseException)
    
end

if ~ischar(controlGainsFilename)
    errID = 'init_sim:BadInput';
    errMsg = 'Expected filename as a string. Filename was %d';
    baseException = MException(errID,errMsg,controlGainsFilename);
    
    causeException = MException(errID, errMsg, controlGainsFilename);
    baseException = addCause(baseException,causeException);
    
    throw(baseException);
end

if exist(controlGainsFilename,'file') ~= 2
    errID = 'init_sim:BadInput';
    errMsg = 'File not found. Filename was %s';
    baseException = MException(errID,errMsg,controlGainsFilename);
    
    causeException = MException(errID, errMsg, controlGainsFilename);
    baseException = addCause(baseException,causeException);
    
    throw(baseException);
end

if exist(iiwaUrdfPath,'file') ~= 2
    errID = 'init_sim:BadInput';
    errMsg = 'File not found. Filename was %s';
    baseException = MException(errID,errMsg,controlGainsFilename);
    
    causeException = MException(errID, errMsg, controlGainsFilename);
    baseException = addCause(baseException,causeException);
    
    throw(baseException);
end

if exist('dt','var') == 0
    errID = 'init_sim:BadInput';
    errMsg = 'dt not found.';
    baseException = MException(errID,errMsg);
    
    causeException = MException(errID, errMsg);
    baseException = addCause(baseException,causeException);
    
    exceptionCorrection = matlab.lang.correction.AppendArgumentsCorrection('0.001');
    baseException = baseException.addCorrection(exceptionCorrection);
    throw(baseException);
end

%%
tmp1 = load(controlGainsFilename);
t_end = length(tmp1.q_mean)*dt;
q_0 = tmp1.q_mean(1,:);

%%
lbr = importrobot(iiwaUrdfPath);
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
end

