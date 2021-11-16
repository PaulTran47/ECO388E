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
% 15 November, 2021
%==========================================================================

%==========================================================================
%% Model info
% y_i = I(U_{i, t} > 0) = I(theta0 + theta1*p_{i, t} + theta2*y_{i, t - 1} + sigma_alpha*\alpha_{i} + e_{i, t} > 0)
% e_{i, t} ~ iid difference between two T1EV deviates (i.e., "logit" error)
% \alpha_i ~ N(0, 1) independently across i
% y_i0 = 0 for all i
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
y = data1(:, 3);
p = data1(:, 4);

% Getting total number of unique consumers.
global N;
N = length(unique(i));

% Getting total time period for each consumer.
global T;
T = length(unique(t));

% Setting total number of simulations used later on when integrating out
% unobserved heterogeneity, \alpha_i. We are choosing 100 draws due to
% parameter estimates being essentially unchanged after testing with 250,
% 500, and 1000 draws.
global S;
S = 100;

clear data1;
cd(home_dir);
%==========================================================================

%==========================================================================
%% Part 3d: Discrete choice model with panel data, state dependence, and unobserved heterogeneity: MLE
%=====
% NOTE
%=====
% Recall we are choosing to integrate out \alpha_i ~ N(0, 1) using a max
% simulated draw number of S = 100. There are N = 100 unique customers.
% The total time period for each customer is T = 20 weeks. We wish for each
% customer to have their own unique simulated draws of unobserved
% heterogeneity. Therefore, \alpha_i is a 3D array of dimension (T, S, N).
% In words, this means each customer has 100 simulated draws of unobserved
% heterogeneity for each time period.

% However, recall that unobserved heterogeneity is time-invariant for each
% consumer i. This means the Tx1 vector for each consumer i and each
% simulated draw s should be the same value.
%=========
% END NOTE
%=========
global alpha;
alpha = randn(1, S, N);
alpha = repmat(alpha, T, 1);

% Reshaping variables y and p into 3D arrays of dimension (T, 1, N). In
% words, each unique customer now has a Tx1 vector of y and p.
global y_it p_it;
y_it = reshape(y, [T, 1, N]);
p_it = reshape(p, [T, 1, N]);

% Creating y_{i, t - 1}. We will also reshape this into a 3D array of
% dimension (T, 1, N). Recall our assumption that y_i0 = 0 for all i.
global y_itm1;
y_itm1 = reshape(lagmatrix(y_it, 1), [T, 1, N]);
y_itm1(isnan(y_itm1)) = 0;

%=====
% NOTE
%=====
% WE ARE CHOOSING TO INTEGRATE OUT \alpha_i.

% Parameters:
% theta(1) = theta0;
% theta(2) = theta1;
% theta(3) = theta2;
% theta(4) = sigma_alpha
%=========
% END NOTE
%=========
% Creating likelihood contribution of i. More information about the
% l_io_alpha_i_i_sum function can be found in l_io_alpha_i_i_sum.m
l_i = @(theta, i) l_io_alpha_i_i_sum(theta, i)/S;
% Creating log-likelihood contribution of i
ll_i = @(theta, i) log(l_i(theta, i));
% Creating sum component of aggregate log-likelihood function. This sums
% across all consumers i.
ll_sum_component = @(theta) 0;
for i = 1:N
  ll_sum_component = @(theta) ll_sum_component(theta) + ll_i(theta, i);
end
clear i;
% Creating aggregate log-likelihood function
ll = @(theta) ll_sum_component(theta)/N;

% Initial parameter values
theta0 = [0.5, 0.5, 0.5, 0.5];

% Numerical optimisation to find (theta1, theta2, sigma)(fminsearch)
%=====
% NOTE
%=====
% After the initial run of the programme to obtain the parameter estimates,
% bhat_ml, we then leave the numerical optimisation code below commented
% out for the rest of the programme's subsequent runs. This is because even
% with S = 100, it still takes a little bit of time to run, which will get
% annoying really quick.
%=========
% END NOTE
%=========
% options = optimset('Display', 'iter');
% bhat_ml = fminsearch(@(theta) -ll(theta), theta0, options);
%=======
% ANSWER
%=======
% We are commenting the estimate of bhat_ml due to the amount of time it
% took for fminsearch to complete with S = 100.

% bhat_ml = [-0.795278688779827, 0.916936802001065, 0.509747530171734, 0.891947440705045].
%         \approx [-0.7953, 0.9169, 0.5097, 0.8919].
%===========
% END ANSWER
%===========

% We are explicitly declaring bhat_ml with the values obtained from
% numerical optimisation. This is so we don't have to run the optimisation
% code again for future parts.
bhat_ml = [-0.795278688779827, 0.916936802001065, 0.509747530171734, 0.891947440705045];

% Calculating the estimated gradients of the log-likelihood contributions
% of i. We are estimating the gradient components numerically by a 
% "one-sided derivative" with steps of 0.001. After testing with
% "two-sided derivatives" in problem 2 of this problem set, the small
% amount of increased accuracy doesn't warrant the increased amount of
% code.

% Making parameter vector to column  vector
bhat_ml = bhat_ml';

%=====
% NOTE
%=====
% Recall that our transformation of the panel data into 3D arrays where
% each z-dimension slice represents a consumer results in the
% (log-)likelihood function to be a scalar for each consumer i. We want our
% gradient components to be vectors, where each row represents a consumer.
%=========
% END NOTE
%=========

% 1st gradient component for i
bhat_ml_theta1 = bhat_ml;
bhat_ml_theta1(1, 1) = bhat_ml(1, 1)*1.001;
dll_idtheta1 = zeros(N, 1);
for i = 1:N
  dll_idtheta1(i, 1) = (ll_i(bhat_ml_theta1, i) - ll_i(bhat_ml, i))/(0.001*bhat_ml(1, 1));
end
clear i bhat_ml_theta1;

% 2nd gradient component for i
bhat_ml_theta2 = bhat_ml;
bhat_ml_theta2(2, 1) = bhat_ml(2, 1)*1.001;
dll_idtheta2 = zeros(N, 1);
for i = 1:N
  dll_idtheta2(i, 1) = (ll_i(bhat_ml_theta2, i) - ll_i(bhat_ml, i))/(0.001*bhat_ml(2, 1));
end
clear i bhat_ml_theta2;

% 3nd gradient component for i
bhat_ml_theta3 = bhat_ml;
bhat_ml_theta3(3, 1) = bhat_ml(3, 1)*1.001;
dll_idtheta3 = zeros(N, 1);
for i = 1:N
  dll_idtheta3(i, 1) = (ll_i(bhat_ml_theta3, i) - ll_i(bhat_ml, i))/(0.001*bhat_ml(3, 1));
end
clear i bhat_ml_theta3;

% 4th gradient component for i
bhat_ml_theta4 = bhat_ml;
bhat_ml_theta4(4, 1) = bhat_ml(4, 1)*1.001;
dll_idtheta4 = zeros(N, 1);
for i = 1:N
  dll_idtheta4(i, 1) = (ll_i(bhat_ml_theta4, i) - ll_i(bhat_ml, i))/(0.001*bhat_ml(4, 1));
end
clear i bhat_ml_theta4;

% Creating numerical gradient for i
dll_idtheta = [dll_idtheta1 dll_idtheta2 dll_idtheta3 dll_idtheta4];
clear dll_idtheta1 dll_idtheta2 dll_idtheta3 dll_idtheta4;
dll_idtheta = dll_idtheta';

% Calculating numerical var-cov matrix
varcov_matrix_ml = zeros(4, 4);
for i = 1:N
  varcov_matrix_ml_i = dll_idtheta(:, i)*dll_idtheta(:, i)';
  varcov_matrix_ml = varcov_matrix_ml + varcov_matrix_ml_i;
end
clear i dll_idtheta varcov_matrix_ml_i;
varcov_matrix_ml = inv(varcov_matrix_ml);
% Obtaining numerical var and SE of (theta0, theta1, theta2, sigma_alpha)
var_bhat_ml = diag(varcov_matrix_ml);
clear varcov_matrix_ml;
se_bhat_ml = sqrt(var_bhat_ml);
clear var_bhat_ml;
%=======
% ANSWER
%=======
% We are commenting the estimate of se_bhat_ml due to the amount of time it
% took for fminsearch to complete with S = 100.

% se_bhat_ml = [0.193604100650191, 0.155842360828226, 0.101160331798504, 0.117196034772911].
%         \approx [0.1936, 0.1558, 0.1012, 0.1172].
%===========
% END ANSWER
%===========
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