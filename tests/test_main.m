
%% Check for ProMP existense
assert(exist('ProMP.m','file') == 2, ...
'ProMP.m file does not exist. did you open the project?')

%% Test other
assert(2 == 2, '2 is not 2')