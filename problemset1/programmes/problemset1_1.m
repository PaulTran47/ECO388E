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

% ECO388E Problem Set 1, 1
% Paul Le Tran, plt377
% 17 September, 2021
%==========================================================================

%==========================================================================
%% Model info
% Model: y_i = theta1 + theta2*x_i + e_i
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
data1 = importdata('data1.dat');

% Setting regressors and regressand
regressors = data1(:, 2);
regressand = data1(:, 1);

% Getting sample size
% Variable is global so it can be accessed in the separate workspaces of
% called functions
global N;
N = length(regressand);

clear data1;
cd(home_dir);
%==========================================================================
%% Part 1a: OLS
% Estimating model
fit = fitlm(regressors, regressand);

% Saving parameters
bhat_ols = fit.Coefficients.Estimate;

% Saving SEs calculated using standard formula
covvar_matrix = fit.CoefficientCovariance;
se_bhat_ols = diag(sqrt(covvar_matrix));
clear fit covvar_matrix;
%==========================================================================

%==========================================================================
%% Part 1b: ML
% Assuming e_i ~ N(0, sigma^2), iid; sigma > 0
% ==> sigma*e_i ~ N(0, 1), iid
% Therefore, model: y_i = theta1 + theta2*x_i + sigma*e_i
% Because there exists an invertible relationship between y_i and e_i, we
% can do ToRV:
% e_i = f^(-1)(y_i, x_i) = y_i - (theta1 - theta2*x_i)/sigma;
% d(f^(-1))/d(y_i) = 1/sigma
% Therefore, p(y_i|x_i; theta1, theta2, sigma) = normcdf(y_i - (theta1 - theta2*x_i)/sigma)*(1/sigma)

% Creating inverse function e_i
% theta(3) = sigma
e_i = @(theta) (regressand - theta(1) - theta(2)*regressors)/abs(theta(3));
% Creating p(y_i|x_i; theta1, theta2, sigma) = likelihood contribution of i
l_i = @(theta) normpdf(e_i(theta))/abs(theta(3));
% Creating log-likelihood contribution of i
ll_i = @(theta) log(l_i(theta));
% Creating log-likelihood function
ll = @(theta) inv(N)*sum(ll_i(theta));

% Initial parameter values for numerical optimisation
theta0 = [1, 1, 1];

% Numerical optimisation to find (theta1, theta2, sigma)(fminsearch)
bhat_ml = fminsearch(@(theta) -ll(theta), theta0);
clear theta0;

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
clear i varcov_matrix_ml_i dll_idtheta;
varcov_matrix_ml = inv(varcov_matrix_ml);

% Obtaining numerical var and SE of (theta1, theta2, sigma)
var_bhat_ml = diag(varcov_matrix_ml);
clear varcov_matrix_ml;
se_bhat_ml = sqrt(var_bhat_ml);
se_bhat_ml = se_bhat_ml(1:2, 1);
clear var_bhat_ml;
%==========================================================================

%==========================================================================
%% Part 1c: GMM
% Redefining inverse function e_i for GMM. Because we have invertibility,
% we don't need to make assumptions on full distribution of e_i. This means
% we don't need to have sigma as a parameter like ML.
e_i = @(theta) regressand - theta(1) - theta(2)*regressors;

% Creating individual moment function g_i (500x2)
g_i = @(theta) [e_i(theta) times(e_i(theta), regressors)];

% Creating mean aggregate moment function GN_i (2x1)
GN = @(theta) inv(N)*sum(g_i(theta))';

% Estimating optimal weight matrix A
% Creating component matrix
% Because this matrix contains the components for every i, it is 2x500
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

% Initial parameter values for numerical optimisation
theta0 = [1, 1];

% Numerical optimisation to find (theta1, theta2)(fminsearch)
bhat_gmm = fminsearch(@(theta) objfunc_gmm(theta), theta0);
clear theta0;

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
%% Questions asked in problem 1 that doesn't involve code
%=========
% Part 1c:
%=========
% Sigma is no longer a parameter that needs to be identified in GMM for
% this model because we have invertibility. This allows us to not have to
% make asssumptions about the full distribution of e_i up to parameters.

% Recall that the asympotatic variance for the GMM estimator is
% inv(G'AG)*(G'AVAG)*inv(G'AG), where G is capital gamma. Furthermore,
% recall that our model is exactly identified since we have both 2 moments
% and 2 parameters. Therefore, we can see that G is 2x2 and V is 2x2.
% From algebra,
% inv(G'AG)*(G'AVAG)*inv(G'AG) = (inv(G)inv(A)inv(G'))*(G'AVAG)*(inv(G)inv(A)inv(G'))
%                              = (inv(G)inv(A))*(AVA)*(inv(A)inv(G'))
% inv(G'AG)*(G'AVAG)*inv(G'AG) = inv(G)(V)inv(G'),
% which is our desired result.
%==========================================================================