
% Known return values from exist(NAME)
RETVAL_EXIST = [];
RETVAL_EXIST.DOES_NOT_EXIST = 0;
RETVAL_EXIST.IS_A_VARIABLE_IN_THE_WORKSPACE = 1;
RETVAL_EXIST.IS_A_FILE = 2;
RETVAL_EXIST.IS_A_MEX_FILE = 3;
RETVAL_EXIST.IS_A_SIMULINK_MODEL_OR_LIB = 4;
RETVAL_EXIST.IS_A_BUILT_IN_MATLAB_FUNCTION = 5;
RETVAL_EXIST.IS_A_P_CODE_FILE = 6;
RETVAL_EXIST.IS_A_FOLDER = 7;
RETVAL_EXIST.IS_A_CLASS = 8;

%% Check for iiwa14_loadMesh.m existence
ERRID = 'MOJO:AssertionFailed:notExist';
ERRMSG = 'iiwa14_loadMesh.m file does not exist. Did you open the project?';
EXPRESSION = exist('iiwa14_loadMesh.m','file') == RETVAL_EXIST.IS_A_FILE;
assert(EXPRESSION, ERRID, ERRMSG)

%% Check for ProMP existence
ERRID = 'MOJO:AssertionFailed:notExist';
ERRMSG = 'ProMP.m file does not exist. Did you open the project?';
EXPRESSION = exist('ProMP.m','file') == RETVAL_EXIST.IS_A_FILE;
assert(EXPRESSION, ERRID, ERRMSG)

%% Test init_sim
ERRID = 'MOJO:AssertionFailed:oopsie';
ERRMSG = 'Oopsie in init_sim.';
EXPRESSION = strcmp(class(init_sim('../controlGains_X-tremeFinal.mat',...
    '../iiwa_description/urdf/iiwa14.urdf', 0.001)),...
    'rigidBodyTree');
assert(EXPRESSION, ERRID, ERRMSG)

%% Test other
ERRID = 'MOJO:AssertionFailed:notExist';
ERRMSG = '2 is not 2';
EXPRESSION = 2 == 2;
assert(EXPRESSION, ERRID, ERRMSG)