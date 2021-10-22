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

% ECO388E Problem Set 1, 2
% Paul Le Tran, plt377
% 19 September, 2021
%==========================================================================

%==========================================================================
%% Model info
% Model: y_i = exp(theta1 + (x_i)^theta2 + sigma*e_i)
%==========================================================================

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'path\to\programmes';
data_dir = 'path\to\data';

cd(home_dir);

%% Loading in data
cd(data_dir);
data2 = importdata('data2.dat');

% Setting regressors and regressand
regressors = data2(:, 2);
regressand = data2(:, 1);

% Getting sample size
% Variable is global so it can be accessed in the separate workspaces of
% called functions
global N;
N = length(regressand);

clear data2;
cd(home_dir);
%==========================================================================

%==========================================================================
%% Part 2a: ML
% Assuming e_i ~ N(0, sigma^2), iid; sigma > 0
% ==> sigma*e_i ~ N(0, 1), iid
% Therefore, model: y_i = exp(theta1 + (x_i)^theta2 + sigma*e_i)
% Because there exists an invertible relationship between y_i and e_i, we
% can do ToRV:
% e_i = (ln(y_i) - theta1 - (x_i)^theta2)/sigma;
% d(f^(-1))/d(y_i) = 1/(y_i*sigma)
% Therefore, p(y_i|x_i; theta1, theta2, sigma) = normpdf((ln(y_i) - theta1 - (x_i)^theta2)/sigma)*(1/(y_i*sigma));

% Creating inverse function e_i
% theta(3) = sigma
e_i = @(theta) (log(regressand) - theta(1) - regressors.^theta(2))/abs(theta(3));
% Creating p(y_i|x_i; theta1, theta2, sigma) = likelihood contribution of i
l_i = @(theta) normpdf(e_i(theta))./abs(regressand*abs(theta(3)));
% Creating log-likelihood contribution of i
ll_i = @(theta) log(l_i(theta));
% Creating log-likelihood function
ll = @(theta) inv(N)*sum(ll_i(theta));

% Because we have a non-linear model, we need to test out different initial
% values to see what the global max is.
% After some manual trial-and-error, I see that any value higher than 8 for
% theta1 and theta2 results in no max to be found. For simplicity, we will
% restrict the domains of the parameters to be
% ([-10, 10], [-10, 10] [1,10]). We will search for the highest local max
% via a loop. We're also skipping by 4 for theta1, theta2 and by 3 for
% theta3/sigma because time.
% Creating counter variable
counter = 0;
% Starting value of the initial parameter values
theta0_start = [-10, -10, 1];
% Creating a matrix that will save every iteration's local max.
bhat_ml_matrix = zeros(144, 3);
for i = theta0_start(1, 1):4:10 % theta1 values
  for j = theta0_start(1, 2):4:10 %theta2 values
    for k = theta0_start(1, 3):3:10 %theta3/sigma values
      theta0 = [i, j, k]; % Actual initial parameter values for numerical optimisation
      % Numerical optimisation to find (theta1, theta2, sigma)(fminsearch)
      bhat_ml = fminsearch(@(theta) -ll(theta), theta0);
      counter = counter + 1;
      bhat_ml_matrix(counter, :) = bhat_ml;
    end
  end
end
clear counter theta0_start i j k theta0 bhat_ml;

% Removing bhat_ml values that are integers equal to the initial values,
% because that indicates no local max was found. Since such situations
% result in all three parameters to be integers, we just need to make logic
% based on the theta1.
condition1 = mod(bhat_ml_matrix(:, 1), 2) == 0;
bhat_ml_matrix(condition1, :) = [];
clear condition1;

% Creating a vector that stores the log-likelihood when given each bhat_ml
ll_vector = zeros(length(bhat_ml_matrix), 1);
for i = 1:length(bhat_ml_matrix)
  ll_vector(i, 1) = ll(bhat_ml_matrix(i, :));
end

% Finding the bhat_ml values which give the maximum log-likelihood.
ll_bhat_ml_matrix = [ll_vector bhat_ml_matrix];
clear ll_vector bhat_ml_matrix;
[m, i] = max(ll_bhat_ml_matrix);
bhat_ml = ll_bhat_ml_matrix(i(1, 1), 2:end);
clear m i ll_bhat_ml_matrix;

% Calculating the estimated gradients of the log-likelihood contributions
% of i. We are estimating the gradient components numerically by a 
% "one-sided derivative" with steps of 0.001.

% Making parameter vector to column  vector
bhat_ml = bhat_ml';

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

% 3nd gradient component for i
bhat_ml_theta3 = bhat_ml;
bhat_ml_theta3(3, 1) = bhat_ml(3, 1)*1.001;
dll_idtheta3 = (ll_i(bhat_ml_theta3) - ll_i(bhat_ml))/(0.001*bhat_ml(3, 1));
clear bhat_ml_theta3;

% Creating numerical gradient for i
dll_idtheta = [dll_idtheta1 dll_idtheta2 dll_idtheta3];
clear dll_idtheta1 dll_idtheta2 dll_idtheta3;
dll_idtheta = dll_idtheta';

% Calculating numerical var-cov matrix
varcov_matrix_ml = zeros(3, 3);
for i = 1:N
  varcov_matrix_ml_i = dll_idtheta(:, i)*dll_idtheta(:, i)';
  varcov_matrix_ml = varcov_matrix_ml + varcov_matrix_ml_i;
end
clear i varcov_matrix_ml_i;
varcov_matrix_ml = inv(varcov_matrix_ml);
% Obtaining numerical var and SE of (theta1, theta2, sigma)
var_bhat_ml = diag(varcov_matrix_ml);
clear varcov_matrix_ml;
se_bhat_ml = sqrt(var_bhat_ml);
se_bhat_ml = se_bhat_ml(1:2, 1);
clear var_bhat_ml;
%==========================================================================

%==========================================================================
%% Part 2b: GMM
% Redefining inverse function e_i for GMM. Because we have invertibility,
% we don't need to make assumptions on full distribution of e_i. This means
% we don't need to have sigma as a parameter like ML.
e_i = @(theta) log(regressand) - theta(1) - regressors.^theta(2);

% Creating individual moment function g_i (100x2)
g_i = @(theta) [e_i(theta) times(e_i(theta), regressors)];

% Creating mean aggregate moment function GN_i (2x1)
GN = @(theta) inv(N)*sum(g_i(theta))';

% Estimating optimal weight matrix A
% Creating component matrix
% Because this matrix contains the components for every i, it is 2x100
% Making function global so its array output can be accessed from all
% workspaces of all called functions.
global component;
component = @(theta) minus(g_i(theta)', GN(theta));

% Creating estimate of optimal weight matrix A
% More information about the function A_component_sum function can be found
% in A_component_sum.m
A = @(theta) inv((inv(N)^2)*A_component_sum(theta));

% Creating GMM objective function
objfunc_gmm = @(theta) GN(theta)'*A(theta)*GN(theta);

% Because we have a non-linear model, we need to test out different initial
% values to see what the global max is.
% For simplicity, we will restrict the domains of the parameters to be
% ([-10, 10], [-10, 10]). We will search for the highest local max
% via a loop. We're also skipping by 4 for theta1, theta2 because time.
% Creating counter variable
counter = 0;
% Starting value of the initial parameter values
theta0_start = [-10, -10];
% Creating a matrix that will save every iteration's local max.
bhat_gmm_matrix = zeros(36, 2);
for i = theta0_start(1, 1):4:10 % theta1 values
  for j = theta0_start(1, 2):4:10 %theta2 values
    theta0 = [i, j]; % Actual initial parameter values for numerical optimisation
    disp(counter);
    disp(theta0);
    % Numerical optimisation to find (theta1, theta2, sigma)(fminsearch)
    bhat_gmm = fminsearch(@(theta) objfunc_gmm(theta), theta0);
    disp(bhat_gmm);
    counter = counter + 1;
    bhat_gmm_matrix(counter, :) = bhat_gmm;
  end
end
clear counter theta0_start i j theta0 bhat_gmm;

% Creating a matrix that stores the moment condition when given each
% bhat_gmm
GN_matrix = zeros(2, length(bhat_gmm_matrix));
for i = 1:length(bhat_gmm_matrix)
  GN_matrix(:, i) = GN(bhat_gmm_matrix(i, :));
end
clear i;
% Though GN is 2x1, we will be transposing the matrix to match with
% bhat_gmm_matrix.
GN_matrix = GN_matrix';

% Finding the bhat_gmm values which give the minimum moment condition.
% To do this, we calculate the Euclidean distance between the moment
% condition vector against the zero vector.
zero_vector = zeros(1, 2);
norm_vector = zeros(36, 1);
for i = 1:length(GN_matrix)
  norm_vector(i, 1) = norm(zero_vector - GN_matrix(i, :));
end
clear zero_vector;
norm_GN_bhat_gmm_matrix = [norm_vector GN_matrix bhat_gmm_matrix];
clear norm_vector GN_matrix bhat_gmm_matrix;
[m, i] = min(norm_GN_bhat_gmm_matrix);
bhat_gmm = norm_GN_bhat_gmm_matrix(i(1, 1), 4:end);
clear m i norm_GN_bhat_gmm_matrix;

% Recall that GN is 2x1 vector. This means if we are taking the derivative
% of GN with respect to bhat_gmm', we are doing vector-by-vector
% derivatives. Our result will be a 2x2 matrix, or the Jacobian matrix.
% The Jacobian matrix is what we call Gamma/G at the bottom of this code.
% Creating Jacobian matrix
J = zeros(2, 2);

% Making bhat_gmm into a column vector
bhat_gmm = bhat_gmm';

% Obtaining the (1, 1) and (2, 1) components
bhat_gmm_theta1 = bhat_gmm;
bhat_gmm_theta1(1, 1) = bhat_gmm(1, 1)*1.001;
dGNdtheta1 = (GN(bhat_gmm_theta1) - GN(bhat_gmm))/(0.001*bhat_gmm(1, 1));
J(1, 1) = dGNdtheta1(1, 1);
J(2, 1) = dGNdtheta1(2, 1);
clear bhat_gmm_theta1 dGNdtheta1;

% Obtaining the (1, 2) and (2, 2) components
bhat_gmm_theta2 = bhat_gmm;
bhat_gmm_theta2(2, 1) = bhat_gmm(2, 1)*1.001;
dGNdtheta2 = (GN(bhat_gmm_theta2) - GN(bhat_gmm))/(0.001*bhat_gmm(2, 1));
J(1, 2) = dGNdtheta2(1, 1);
J(2, 2) = dGNdtheta2(2, 1);
clear bhat_gmm_theta2 dGNdtheta2;
Gamma = abs(J);

% Obtaining GMM estimator var-cov matrix
varcov_matrix_gmm = inv(Gamma)*inv(A(bhat_gmm))*inv(Gamma');
clear J Gamma V;
% Obtaining numerical var and SE of (theta1, theta2)
var_bhat_gmm = diag(varcov_matrix_gmm);
clear varcov_matrix_gmm;
se_bhat_gmm = sqrt(var_bhat_gmm);
clear var_bhat_gmm;
%==========================================================================

%==========================================================================
%% Questions asked in problem 2 that doesn't involve code
%=======================
% (Technically) Part 2c:
%=======================
% Recall that for ML, you need to assume the full distribution of the error
% term up to parameters. This includes making a stance on what the variance
% of the error terms are. In contrast, for invertible GMM, you only need to
% assume the mean independence, E[e_i|x_i] = 0. This means that invertible
% GMM isn't making a stance on the variance of the error term, and is
% therefore robust to heteroscedasticity.

%=======================
% (Technically) Part 2d:
%=======================
% In the context of this problem, recall that the model is exactly
% identified because the number of moments is equal to the number of
% parameters estimated. As a result, the choice of weight matrix A should
% not matter. We know that when it comes to non-invertible GMM estimation,
% the estimator is not efficient because it does not contain all possible
% higher-order moments of \epsilon_i’s distribution, which capture more
% information about the variance). This is even the case when we make an
% assumption on the distribution of \epsilon_i up to parameters.

%When considering invertible GMM estimation, we only make assume
% E[e_i | x_i] = 0, meaning we do not know anything else about e_i's
% distribution. Therefore, even with optimal weight matrix A and the model
% exactly identified, we are unable to add any high-order moments to the
% moment condition due to not knowing anything else about e_i’s
% distribution. It is this lack of higher-order moments that makes both
% invertible GMM estimators (this problem) and non-invertible GMM
% estimators (problem 3) less efficient than ML estimators.
%==========================================================================