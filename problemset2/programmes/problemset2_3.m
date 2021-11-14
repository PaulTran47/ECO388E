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

% ECO388E Problem Set 2, 3
% Paul Le Tran, plt377
% 14 November, 2021
%==========================================================================

%==========================================================================
%% Model info
% y_i = I(U_{i, t} > 0) = I(theta0 + theta1*p_{i, t} + theta2*y_{i, t - 1} + sigma_alpha*\alpha_{i} + e_{i, t} > 0)
% e_{i, t} ~ iid difference between two T1EV deviates (i.e., "logit" error)
% \alpha_i ~ N(0, 1) independently across i.
%==========================================================================

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'path\to\programmes';
data_dir = 'path\to\data';

cd(home_dir);

%% Loading in simulated data on a panel of consumer decisions whether or not to purchase a grocery item in a given week.
cd(data_dir);
data1 = importdata('sim3.dat');

% Declaring variables
i = data1(:, 1);
t = data1(:, 2);
y_it = data1(:, 3);
p_it = data1(:, 4);

% Getting total number of consumers. This is found by looking at the vector
% i
global N;
N = 100;

% Getting total time period for each consumer. This is found by looking at
% the vector t
global T;
T = 20;

clear data1;
cd(home_dir);
%==========================================================================

%==========================================================================
%% Questions asked in problemset 2_3 that don't involve code
%=========
% Part 3a:
%=========
% Observe that our model has both state dependence and unobserved
% heterogeneity. Because there is correlation in each individual's
% decisions over time, MLE needs to use the joint probability of consumer
% i's decisions in all time periods (conditioned on prices,
% initial decision, and parameters). This means the likelihood of the data
% for consumer i is:

% L_i = p(y_{i, 1}, \ldots, y_{i, T} | p_{i, 1}, \ldots, p_{i, T}, y_{i, 0}; theta).

% Integrating out \alpha_i, our likelihood becomes:

% L_i = \int p(y_{i, 1}, \ldots, y_{i, T} | p_{i, 1}, \ldots, p_{i, T}, y_{i, 0}, \alpha_i; theta)*p(\alpha_i).

% Observe that U_{i, t} only depends on e_{i, t} if we know p_{i, t},
% y_{i, t - 1}, and \alpha_{i}. Furthermore, focus on the probability

% p(I(U_{i, t} > 0) | p_{i, t}, y_{i, t - 1}, \alpha_i; theta).

% Because e_{i, t} are iid, this means they are independent over time.
% Observe that this holds because the individual unobserved heterogeneity
% is integrated out nad accounted for.
% As a result, we know the following holds:

% p(I(U_{i, t} > 0) | p_{i, t}, y_{i, t - 1}, \alpha_i; theta) \perp p(I(U_{i, t + 1} > 0) | p_{i, t + 1}, y_{i, t}, \alpha_i; theta).

% Consider we also know that y_{i, t} = I(U_{i, t} > 0), our relation
% simplifies as follows:

% p(y_{i, t} | p_{i, t}, y_{i, t - 1}, \alpha_i; theta) \perp p(y_{i, t + 1} | p_{i, t + 1}, y_{i, t}, \alpha_i; theta).

% In other words, because e_{i, t} is iid, when combined with conditioning
% on p_{i, t}, y_{i, t - 1}, \alpha_i, and parameters, factoring the inner
% joint probability of our likelihood function for consumer i becomes

% L_i = \int \[\prod_{t = 1}^{T} p(y_{i, t} | p_{i, t}, y_{i, t - 1}, \alpha_i; theta)\]*p(\alpha_i).

%=========
% Part 3b:
%=========
% Consider the following likelihood function:

% L_i = p(y_{i, 1}, \ldots, y_{i, T} | p_{i, 1}, \ldots, p_{i, T}, y_{i, 0}; theta) = \prod_{t = 1}^{T} p(y_{i, t} | p_{i, t}, y_{i, t - 1}; theta).

% This equality would not hold, meaning it would not be obtainable if we
% choose not to integrate out \alpha_i.

% The reason for this is because if we choose not to integrate out the
% unobserved heterogeneity, the e_{i, t} over time could be correlated with
% each other. This would result in the probabilities inside of the product
% to no longer be independent of each other.

%=========
% Part 3c:
%=========
% Recall that y_i = I(U_{i, t} > 0). That means the inner probability terms
% of the product is:

% p(y_{i, t} | p_{i, t}, y_{i, t - 1}, \alpha_i; theta) = p(I(U_{i, t} > 0) | p_{i, t}, y_{i, t - 1}, \alpha_i; theta).

% When we observe y_{i, t} = 1, we know U_{i, t} > 0. This means we can
% write the RHS probability above as:

% p(I(U_{i, t} > 0) | p_{i, t}, y_{i, t - 1}, \alpha_i; theta) = \frac{exp(theta0 + theta1*p_{i, t} + theta2*y_{i, t - 1} + sigma_alpha*\alpha_{i})}{1 + exp(theta0 + theta1*p_{i, t} + theta2*y_{i, t - 1} + sigma_alpha*\alpha_{i})}.

% Similarly, when we observe y_{i, t} = 0, we know U_{i, t} < 0. Our RHS
% probability then becomes:

% p(I(U_{i, t} < 0) | p_{i, t}, y_{i, t - 1}, \alpha_i; theta) = \frac{1}{1 + exp(theta0 + theta1*p_{i, t} + theta2*y_{i, t - 1} + sigma_alpha*\alpha_{i})}.

% Combining the two, the inner probability of the product in our likelihood
% function can therefore be expressed as:

% p(y_{i, t} | p_{i, t}, y_{i, t - 1}, \alpha_i; theta) = y_{i, t}*\frac{exp(theta0 + theta1*p_{i, t} + theta2*y_{i, t - 1} + sigma_alpha*\alpha_{i})}{1 + exp(theta0 + theta1*p_{i, t} + theta2*y_{i, t - 1} + sigma_alpha*\alpha_{i})} + (1 - y_{i, t})*\frac{1}{1 + exp(theta0 + theta1*p_{i, t} + theta2*y_{i, t - 1} + sigma_alpha*\alpha_{i})}.

% Substituting this into our inner-probabability-product will give us
% expression (3) of the likelihood function displayed in the problem set.
%==========================================================================