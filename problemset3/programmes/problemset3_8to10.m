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

% ECO388E Problem Set 3; 8, 9, 10
% Paul Le Tran, plt377
% 12 December, 2021
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

% Creating array that has a and i as columns in that order
data = [a i];

% Creating variable that houses total number of observations
global N;
N = length(a);
%==========================================================================

%==========================================================================
%% Part 8: Rust (1987) - Pseudo-likelihood estimation for parameters (\theta_{1}, R) using Hotz-Miller (1993)
cd(home_dir);
%=====
% NOTE
%=====
% The Hotz-Miller (1993) (HM) approach uses the following alternative
% representative equations for the value function and policy function:

% V(.; \theta) = h(P(.; \theta); \theta);
% P(.; \theta) = g(V(.; \theta); \theta).

% Substituting the first into the second yields the following expression
% for the policy function:

% P(.; \theta) = g(h(P(.; \theta); \theta); \theta).

% The HM approach then replaces P(.) on the RHS of the equation with
% consistent estimator, \hat{P(.)}. This will be estimated
% non-parametrically by simply computing the probabilities of various
% actions found in the data.
%=========
% END NOTE
%=========
% Initialising discounting parameter and Euler's constant
global beta gamma;
beta = 0.9;
gamma = 0.5772;

% Creating variable that houses the maximum age a machine can be (recall a
% machine can only have an age in the set (1:5))
global a_max;
a_max = 5;

% Calculating observed probability of replacing machine conditional on age
p1_hat = splitapply(@sum, i, a)./groupcounts(a);

% Calculating observed probability of keeping machine conditional on age
p0_hat = 1 - p1_hat;

% Creating matrices of transition probabilities, T0 and T1. These represent
% p(a_{t + 1} | a_{t}, i_{t} = 0; \theta) and
% p(a_{t + 1} | a_{t}, i_{t} = 1; \theta) respectively
T1 = zeros(a_max);
T1(:, 1) = ones(a_max, 1);
T0 = vertcat(eye(a_max - 1), zeros(1, a_max - 1));
T0 = horzcat(zeros(a_max, 1), T0);
T0(a_max, a_max) = 1;

%=================================================
% Calculating value function using the HM approach
%=================================================
%=====
% NOTE
%=====
% Paramater value mapping:
% theta(1) = \theta_{1}
% theta(2) = R;
%=========
% END NOTE
%=========
% Creating vector of possible ages
a_hm = (1:5)';

% Creating an a_max x a_max matrix that repeats column of p1_hat
p1_hat_matrix = repmat(p1_hat, 1, a_max);

% Creating an a_max x a_max matrix that repeats column of p0_hat
p0_hat_matrix = repmat(p0_hat, 1, a_max);

% Creating inverse matrix component of value function using the HM approach
inv_component = inv(eye(a_max) - beta.*(p0_hat_matrix).*T0 - beta.*p1_hat_matrix.*T1);

% Creating function that is the non-inverse matrix component of the value
% function using the HM approach
noninv_component = @(theta) (p0_hat).*(theta(1).*a_hm + gamma - log(p0_hat)) + p1_hat.*(theta(2) + gamma - log(p1_hat));

% Creating function for value function using the HM approach
v_hm = @(theta) inv_component*noninv_component(theta);

%==================================================
% Calculating policy function using the HM approach
%==================================================
%=====
% NOTE
%=====
% Paramater value mapping:
% theta(1) = \theta_{1}
% theta(2) = R;
%=========
% END NOTE
%=========
% Using the class notes notation, we will call this function/expression as
% \kappa

% Creating function for numerator of \kappa
numer = @(theta) exp(theta(2) + beta.*T1*v_hm(theta));

% Creating function for denominator of \kappa
denom = @(theta) exp(theta(1).*a_hm + beta.*T0*v_hm(theta)) + numer(theta);

% Creating \kappa function
kappa = @(theta) numer(theta)./denom(theta);

%==================================================================
% Pseudo-likelihood estimation for parameters (\theta_{1}, R) using
% Hotz-Miller (1993)
%==================================================================
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

% Creating pseudo-log-likelihood function for firm i
ll_pseudo_i = @(theta) (1 - i).*log(1 - subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}}))) + i.*log(subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}})));

% Creating objective function
ll_pseudo = @(theta) sum(ll_pseudo_i(theta))/N;

% Numerical optimisation to find (theta1, R)(fminsearch)
% options = optimset('Display', 'iter');
bhat_ml_pseudo = fminsearch(@(theta) -ll_pseudo(theta), theta0);
clear theta0;

% Saving the parameter estimates from problem 5
bhat_ml = [-1.14838831774149, -4.44640773763608];

% Creating vector that will store difference between parameter estimates
% from problem 5 and problem 8. This will also be used to calculate the
% difference between parameter estimates from problem 5 and the upcoming
% problem 9.
error_5_89 = zeros(3, 2);
error_5_89(1, :) = bhat_ml - bhat_ml_pseudo;
%=======
% ANSWER
%=======
% Using pseudo-log-likelihood to estimate our parameters (\theta_{1}, R),
% we obtain the following values for (\hat{\theta_{1}}, \hat{R}) =
% (-1.14950668226087, -4.45374809397243). Recall that our estimates using
% value function iteration are (-1.14838831774149, -4.44640773763608). We
% can see that the two estimates are quite close to each other already, and
% only very slightly off.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 9: Making estimates from HM approach converge to estimates using value function iteration and MLE via the AM iteration process
%=====
% NOTE
%=====
% Due to estiamtes from problem 8 already being quite close to those
% obtained using MLE and value function iteration in problem 5, we will
% iterate the HM approach manually. This is because we assume that it will
% take very little number of iterations for convergence.

% The iteration process itself involves using the HM approach to obtain
% parameter estimates, updating replacing observed p1_hat and p0_hat with
% probabilities obtained from \kappa function when plugging in parameter
% estimates, then repeating the process.
%=========
% END NOTE
%=========
%============
% Iteration 1
%============
% Updating p1_hat, p0_hat, and associated matrices
p1_hat = kappa(bhat_ml_pseudo);
p0_hat = 1 - p1_hat;
p1_hat_matrix = repmat(p1_hat, 1, a_max);
p0_hat_matrix = repmat(p0_hat, 1, a_max);

% Updating all functions
inv_component = inv(eye(a_max) - beta.*(p0_hat_matrix).*T0 - beta.*p1_hat_matrix.*T1);
noninv_component = @(theta) (p0_hat).*(theta(1).*a_hm + gamma - log(p0_hat)) + p1_hat.*(theta(2) + gamma - log(p1_hat));
v_hm = @(theta) inv_component*noninv_component(theta);
numer = @(theta) exp(theta(2) + beta.*T1*v_hm(theta));
denom = @(theta) exp(theta(1).*a_hm + beta.*T0*v_hm(theta)) + numer(theta);
kappa = @(theta) numer(theta)./denom(theta);
ll_pseudo_i = @(theta) (1 - i).*log(1 - subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}}))) + i.*log(subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}})));
ll_pseudo = @(theta) sum(ll_pseudo_i(theta))/N;

% Numerical optimisation to find (theta1, R)(fminsearch)
% options = optimset('Display', 'iter');
bhat_ml_pseudo_1 = fminsearch(@(theta) -ll_pseudo(theta), bhat_ml_pseudo);

% Calculating error
error_5_89(2, :) = bhat_ml - bhat_ml_pseudo_1;

%============
% Iteration 2
%============
% Updating p1_hat, p0_hat, and associated matrices
p1_hat = kappa(bhat_ml_pseudo_1);
p0_hat = 1 - p1_hat;
p1_hat_matrix = repmat(p1_hat, 1, a_max);
p0_hat_matrix = repmat(p0_hat, 1, a_max);

% Updating all functions
inv_component = inv(eye(a_max) - beta.*(p0_hat_matrix).*T0 - beta.*p1_hat_matrix.*T1);
noninv_component = @(theta) (p0_hat).*(theta(1).*a_hm + gamma - log(p0_hat)) + p1_hat.*(theta(2) + gamma - log(p1_hat));
v_hm = @(theta) inv_component*noninv_component(theta);
numer = @(theta) exp(theta(2) + beta.*T1*v_hm(theta));
denom = @(theta) exp(theta(1).*a_hm + beta.*T0*v_hm(theta)) + numer(theta);
kappa = @(theta) numer(theta)./denom(theta);
ll_pseudo_i = @(theta) (1 - i).*log(1 - subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}}))) + i.*log(subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}})));
ll_pseudo = @(theta) sum(ll_pseudo_i(theta))/N;

% Numerical optimisation to find (theta1, R)(fminsearch)
% options = optimset('Display', 'iter');
bhat_ml_pseudo_2 = fminsearch(@(theta) -ll_pseudo(theta), bhat_ml_pseudo_1);

% Calculating error
error_5_89(3, :) = bhat_ml - bhat_ml_pseudo_2;
%=======
% ANSWER
%=======
% After two iterations, we see that our parameter estimates for
% (\theta_{1}, R), (\hat{\theta_{1}}, \hat{R}), are
% (-1.14839273001228, -4.44646620719735). This is nearly identical to the
% parameter estimates obtained by MLE and value function iteration in
% problem 5.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 10: Making estimates from HM approach converge to estimates using value function iteration and MLE via the AM iteration process and random initial probabilities
%============
% Iteration 0
%============
% Making p1_hat begin off as random probabilities drawn from U(0, 1)
p1_hat = rand(a_max, 1);

% Updating p0_hat and associated matrices
p0_hat = 1 - p1_hat;
p1_hat_matrix = repmat(p1_hat, 1, a_max);
p0_hat_matrix = repmat(p0_hat, 1, a_max);

% Updating all functions
inv_component = inv(eye(a_max) - beta.*(p0_hat_matrix).*T0 - beta.*p1_hat_matrix.*T1);
noninv_component = @(theta) (p0_hat).*(theta(1).*a_hm + gamma - log(p0_hat)) + p1_hat.*(theta(2) + gamma - log(p1_hat));
v_hm = @(theta) inv_component*noninv_component(theta);
numer = @(theta) exp(theta(2) + beta.*T1*v_hm(theta));
denom = @(theta) exp(theta(1).*a_hm + beta.*T0*v_hm(theta)) + numer(theta);
kappa = @(theta) numer(theta)./denom(theta);
ll_pseudo_i = @(theta) (1 - i).*log(1 - subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}}))) + i.*log(subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}})));
ll_pseudo = @(theta) sum(ll_pseudo_i(theta))/N;

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

% Numerical optimisation to find (theta1, R)(fminsearch)
% options = optimset('Display', 'iter');
bhat_ml_pseudo_p1_hat_rand = fminsearch(@(theta) -ll_pseudo(theta), theta0);
clear theta0;

% Creating vector that will store difference between parameter estimates
% from problem 5 and problem 10.
error_5_10 = zeros(3, 2);
error_5_10(1, :) = bhat_ml - bhat_ml_pseudo_p1_hat_rand;

%============
% Iteration 1
%============
% Updating p1_hat, p0_hat, and associated matrices
p1_hat = kappa(bhat_ml_pseudo_p1_hat_rand);
p0_hat = 1 - p1_hat;
p1_hat_matrix = repmat(p1_hat, 1, a_max);
p0_hat_matrix = repmat(p0_hat, 1, a_max);

% Updating all functions
inv_component = inv(eye(a_max) - beta.*(p0_hat_matrix).*T0 - beta.*p1_hat_matrix.*T1);
noninv_component = @(theta) (p0_hat).*(theta(1).*a_hm + gamma - log(p0_hat)) + p1_hat.*(theta(2) + gamma - log(p1_hat));
v_hm = @(theta) inv_component*noninv_component(theta);
numer = @(theta) exp(theta(2) + beta.*T1*v_hm(theta));
denom = @(theta) exp(theta(1).*a_hm + beta.*T0*v_hm(theta)) + numer(theta);
kappa = @(theta) numer(theta)./denom(theta);
ll_pseudo_i = @(theta) (1 - i).*log(1 - subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}}))) + i.*log(subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}})));
ll_pseudo = @(theta) sum(ll_pseudo_i(theta))/N;

% Numerical optimisation to find (theta1, R)(fminsearch)
% options = optimset('Display', 'iter');
bhat_ml_pseudo_p1_hat_rand_1 = fminsearch(@(theta) -ll_pseudo(theta), bhat_ml_pseudo_p1_hat_rand);

% Calculating error
error_5_10(2, :) = bhat_ml - bhat_ml_pseudo_p1_hat_rand_1;

%============
% Iteration 2
%============
% Updating p1_hat, p0_hat, and associated matrices
p1_hat = kappa(bhat_ml_pseudo_p1_hat_rand_1);
p0_hat = 1 - p1_hat;
p1_hat_matrix = repmat(p1_hat, 1, a_max);
p0_hat_matrix = repmat(p0_hat, 1, a_max);

% Updating all functions
inv_component = inv(eye(a_max) - beta.*(p0_hat_matrix).*T0 - beta.*p1_hat_matrix.*T1);
noninv_component = @(theta) (p0_hat).*(theta(1).*a_hm + gamma - log(p0_hat)) + p1_hat.*(theta(2) + gamma - log(p1_hat));
v_hm = @(theta) inv_component*noninv_component(theta);
numer = @(theta) exp(theta(2) + beta.*T1*v_hm(theta));
denom = @(theta) exp(theta(1).*a_hm + beta.*T0*v_hm(theta)) + numer(theta);
kappa = @(theta) numer(theta)./denom(theta);
ll_pseudo_i = @(theta) (1 - i).*log(1 - subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}}))) + i.*log(subsref(kappa(theta), struct('type', '()', 'subs', {{a, 1}})));
ll_pseudo = @(theta) sum(ll_pseudo_i(theta))/N;

% Numerical optimisation to find (theta1, R)(fminsearch)
% options = optimset('Display', 'iter');
bhat_ml_pseudo_p1_hat_rand_2 = fminsearch(@(theta) -ll_pseudo(theta), bhat_ml_pseudo_p1_hat_rand_1);

% Calculating error
error_5_10(3, :) = bhat_ml - bhat_ml_pseudo_p1_hat_rand_2;

%=======
% ANSWER
%=======
% The key difference between this problem and problem 9 is that our initial
% value for conditional probabilities p1_hat are not from observed data,
% but instead chosen randomly. Despite this difference, we find that still
% after two iterations, our parameter estimates for
% (\theta_{1}, R), (\hat{\theta_{1}}, \hat{R}), are
% (-1.14838916350878, -4.44642925238419). This is nearly identical to the
% parameter estimates obtained by MLE and value function iteration in
% problem 5.
%===========
% END ANSWER
%===========
%==========================================================================