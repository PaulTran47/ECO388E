%==========================================================================
%% 2 percentage signs represent sections of code;
% 1 percentage sign represents comments for code or commented out code;

% Answers to question parts that don't involve code can be found at the
% bottom of the programme, in the section ``Questions asked in problemset x
% that don't involve code".

% Text answers to question parts that involve code will be between the
% sub-section label:
%=======
% ANSWER
%=======
% Answer here
%===========
% END ANSWER
%===========

% Comments that are important will be between the sub-section label:
%=====
% NOTE
%=====
% Important note here
%=========
% END NOTE
%=========

% ECO388E Problem Set 3; 5, 6
% Paul Le Tran, plt377
% 11 December, 2021
%==========================================================================

%==========================================================================
%% Model info
% Similar capital replacement value function iteration problem as in Rust
% (1987)
% Firms each use one machine to produce output in each period. These
% machines age and become more likely to break down over time. Each period,
% firms have option to replace the machines.
% a_{t} = age of machine at time t
% \Pi(a_{t}, i_{t}, \epsilon_{0, t}, \epsilon_{1, t}) = 
%   \theta_{1}a_{t} + \epsilon_{0, t}, if i_{t} = 0;
%   R + \epsilon_{1, t}, if i_{t} = 1,
% where i_{t} = 1 if firm decides to replace machine at t, R = net cost of
% new machine, and \epsilon_{t}'s are time specific shocks to the utilities
% from replacing and not replacing. Assume that \epsilon_{t}'s are iid
% logit errors.
% State evolution equation:
% a_{t + 1} = 
%  min(5, a_{t} + 1), if i_{t} = 0;
%  1, if i_{t} = 1.
% Note: After 5 years machines don't age. a_{t} \in {1, 2, 3, 4, 5}.
%==========================================================================

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'path\to\programmes';
data_dir = 'path\to\data';

cd(home_dir);

% Loading in data
cd(data_dir);
filename = append(data_dir, '\data.asc');
% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: categorical (%C)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%16f%17f%C%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
% Close the text file.
fclose(fileID);
% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post processing code is included. To generate code which works for unimportable data, select unimportable cells in a file and regenerate the script.
% Create output variable
data = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','VarName2','VarName3'});
% Clear temporary variables
clearvars filename formatSpec fileID dataArray ans;

% Extracting columns from table to make vectors. In terms of meaning,
% consider this as cross-sectional data. This means that each row is a
% different firm.
a = data{:, 1};
i = data{:, 2};
clear data;

% Creating variable that houses total number of observations
N = length(a);
%==========================================================================

%==========================================================================
%% Part 5: Rust (1987) - MLE for parameters (\theta_{1}, R)
cd(home_dir);

% Initialising discounting parameter and Euler's constant
global beta gamma;
beta = 0.9;
gamma = 0.5772;

% Creating variable that houses the maximum age a machine can be (recall a
% machine can only have an age in the set (1:5))
global a_max;
a_max = 5;

% Creating initial value for parameters
%=====
% NOTE
%=====
% Paramater value mapping:
% theta(1) = \theta_{1}
% theta(2) = R;
%=========
% END NOTE
%=========
theta0 = [-1 -3];

% Performing value function iteration. More information about the vfi
% function can be found in vfi.m
v0v1_matrix = @(theta) vfi(theta);

% Creating functions for p(i = 1 | a_{t}) and p(i = 0 | a_{t})
p0 = @(theta) exp(subsref(v0v1_matrix(theta), struct('type', '()', 'subs', {{a, 1}})))./(exp(subsref(v0v1_matrix(theta), struct('type', '()', 'subs', {{a, 1}}))) + exp(subsref(v0v1_matrix(theta), struct('type', '()', 'subs', {{a, 2}}))));
p1 = @(theta) 1 - p0(theta);

% Creating likelihood function for firm i
l_i = @(theta) (1 - i).*p0(theta) + i.*p1(theta);

% Creating log-likelihood function for firm i
ll_i = @(theta) log(l_i(theta));

% Creating objective function
ll = @(theta) sum(ll_i(theta))/N;

% Numerical optimisation to find (theta1, R)(fminsearch)
% options = optimset('Display', 'iter');
bhat_ml = fminsearch(@(theta) -ll(theta), theta0);
clear theta0;

%=======
% ANSWER
%=======
% Our MLE estimation for parameters (\theta_{1}, R),
% (\hat{\theta_{1}}, \hat{R}), is (-1.14838831774149, -4.44640773763608).

% Additionally, we are estimating based on the likelihood of i_{t} being
% conditional on a_{t} because the firm's decision to replace the machine
% depends on how old the machine is. Specifically, the firm's decision is
% based on the expected future costs of its machine. This in turn is
% dependent on the current age of the machine. As a result, the likelihood
% must be condition on age a_{t} in order to do our estimation of
% parameters.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 6: Rust (1987) - Calculating SE for (\hat{\theta_{1}}, \hat{R})
% Calculating the estimated gradients of the log-likelihood contributions
% of i. We are estimating the gradient components numerically by a 
% "one-sided derivative" with steps of 0.001. After testing with
% "two-sided derivatives" in problem 2 of this problem set, the small
% amount of increased accuracy doesn't warrant the increased amount of
% code.

% Making parameter vector to column vector
bhat_ml = bhat_ml';

%=====
% NOTE
%=====
% Recall that our transformation of the panel data into 3D arrays where
% each z-dimension slice represents a firm results in the
% (log-)likelihood function to be a scalar for each firm i. We want our
% gradient components to be vectors, where each row represents a firm.
%=========
% END NOTE
%=========

% 1st gradient component for i
bhat_ml_theta1 = bhat_ml;
bhat_ml_theta1(1, 1) = bhat_ml(1, 1)*1.001;
dll_idtheta1 = (ll_i(bhat_ml_theta1) - ll_i(bhat_ml))/(0.001*bhat_ml(1, 1));
clear bhat_ml_theta1;

% 2nd gradient component for i
bhat_ml_theta2 = bhat_ml;
bhat_ml_theta2(2, 1) = bhat_ml(2, 1)*1.001;
dll_idtheta2 = (ll_i(bhat_ml_theta2) - ll_i(bhat_ml))/(0.001*bhat_ml(2, 1));
clear bhat_ml_theta2;

% Creating numerical gradient for i
dll_idtheta = [dll_idtheta1 dll_idtheta2];
clear dll_idtheta1 dll_idtheta2;
dll_idtheta = dll_idtheta';

% Calculating numerical var-cov matrix
varcov_matrix_ml = zeros(2, 2);
for i = 1:N
  varcov_matrix_ml_i = dll_idtheta(:, i)*dll_idtheta(:, i)';
  varcov_matrix_ml = varcov_matrix_ml + varcov_matrix_ml_i;
end
clear i varcov_matrix_ml_i dll_idtheta;
varcov_matrix_ml = inv(varcov_matrix_ml);

% Obtaining numerical var and SE of (theta1, R)
var_bhat_ml = diag(varcov_matrix_ml);
clear varcov_matrix_ml;
se_bhat_ml = sqrt(var_bhat_ml);

%=======
% ANSWER
%=======
% The SEs for our MLE of parameters (\theta_{1}, R),
% ((SE(\hat{\theta_{1}}), \hat{R})), is
% (0.0748766496059878, 0.320982104639120).
%===========
% END ANSWER
%===========
%==========================================================================