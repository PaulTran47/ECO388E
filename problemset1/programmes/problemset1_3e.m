%==========================================================================
%% 2 percentage signs represent sections of code;
% 1 percentage sign represents comments for code or commented out code;

% ECO388E Problem Set 1, 3e
% Paul Le Tran, plt377
% 21 September, 2021
%==========================================================================

%==========================================================================
%% Model info (3a - 3e)
% q_i = q_1 + q_2
% p_i = 100 - q_i = 100 - q_1 - q_2
% mc_1 = exp(theta1 + theta2*x1_i + sigma*e1_i)
% mc_2 = exp(theta1 + theta2*x2_i + sigma*e2_i)
% From 3a, we get the following model:
% p_i = (100 + exp(theta1 + theta2*x1_i + sigma*e1_i) + exp(theta1 + theta2*x2_i + sigma*e2_i))/3
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
data3 = importdata('data3.dat');
data_drawsgmm = importdata('drawsgmm.dat');

% Setting regressors and regressand
% Variables are global so they can be accessed in the separate workspaces
% of called functions
global x1_i x2_i regressand;
x1_i = data3(:, 2);
x2_i = data3(:, 3);
regressand = data3(:, 1);

% Getting sample size (in this case, sample size represents number of
% markets both firms are in)
% Variable is global so it can be accessed in the separate workspaces of
% called functions
global N;
N = length(regressand);

clear data3;
cd(home_dir);
%==========================================================================

%==========================================================================
%% Part 3e: GMM (Technically MSM)
% Observe that by removing the (p_i)^2 moment, our model is now exactly
% specified. This is because we have 3 moments and 3 parameters.

% Creating matrix that houses all the simulated draws for e1_i and e2_i
% Both have 20 draws, so the matrix is 50x40 due to N = 50. The first 20
% columns are draws for e1_i. The second 20 draws are draws for e2_i.
% Matrix is global so it can be accessed in the separate workspaces of
% called functions
global e_i_matrix;
e_i_matrix = data_drawsgmm;
clear data_drawsgmm;

% Creating variable storing number of simulated draws
% Variable is global so it can be accessed in the separate workspaces of
% called functions
global S;
S = min(size(e_i_matrix))/2;

% Creating the left component of the first moment in the individual moment
% function g_i. This is the difference between p_i and the simulated inner
% expectations of p_i (50x1). More info about the inner_expectations_p_sum
% function can be found in inner_expectations_p_sum.m
g_1st_m_component_i = @(theta) regressand - inv(S)*inner_expectations_p_sum(theta);

% Creating the individual moment function g_i (50x3)
g_i = @(theta) [g_1st_m_component_i(theta) g_1st_m_component_i(theta).*x1_i g_1st_m_component_i(theta).*x2_i];

% Creating mean aggregate moment function GN_i (3x1)
GN = @(theta) inv(N)*sum(g_i(theta))';

% Estimating optimal weight matrix A
% Creating component matrix
% Because this matrix contains the components for every i, it is 4x50
% Making function global so its array output can be accessed from all
% workspaces of all called functions.
global component;
component = @(theta) minus(g_i(theta)', GN(theta));

% Creating estimate of optimal weight matrix A
% More information about the function A_component_sum function can be found
% in A_component_sum.m
A = @(theta) inv((inv(N)^2)*A_component_sum(theta));

%=============================
% FIRST STAGE ESTIMATION START
%=============================

% Creating GMM objective function for first stage
objfunc_gmm_fs = @(theta) GN(theta)'*A(theta)*GN(theta);

% For first stage, we will use the initial parameter values of [1, 1, .01].
% We will then save the parameter values that are spit out from the
% numerical optimisation.
theta0 = [1, 1, 0.01];
bhat_gmm_fs = fminsearch(@(theta) objfunc_gmm_fs(theta), theta0);
clear theta0;

% Calculating the A matrix using bhat_gmm_fs found in the first stage.
% This will be used for creating the first stage varcov matrix and in
% second stage estimation.
A_bhat_gmm_fs = A(bhat_gmm_fs);

% Recall that GN is 3x1 vector. This means if we are taking the derivative
% of GN with respect to bhat_gmm', we are doing vector-by-vector
% derivatives. Our result will be a 3x3 matrix, or the Jacobian matrix.
% The Jacobian matrix is what we call Gamma/G at the bottom of this code.
% Creating Jacobian matrix
J_fs = zeros(3, 3);

% Making bhat_gmm_fs into a column vector
bhat_gmm_fs = bhat_gmm_fs';

% Obtaining the (1, 1), (2, 1), and (3, 1) components
bhat_gmm_fs_theta1 = bhat_gmm_fs;
bhat_gmm_fs_theta1(1, 1) = bhat_gmm_fs(1, 1)*1.001;
dGNdtheta1_fs = (GN(bhat_gmm_fs_theta1) - GN(bhat_gmm_fs))/(0.001*bhat_gmm_fs(1, 1));
J_fs(1, 1) = dGNdtheta1_fs(1, 1);
J_fs(2, 1) = dGNdtheta1_fs(2, 1);
J_fs(3, 1) = dGNdtheta1_fs(3, 1);
clear bhat_gmm_fs_theta1 dGNdtheta1_fs;

% Obtaining the (1, 2), (2, 2), and (3, 2) components
bhat_gmm_fs_theta2 = bhat_gmm_fs;
bhat_gmm_fs_theta2(2, 1) = bhat_gmm_fs(2, 1)*1.001;
dGNdtheta2_fs = (GN(bhat_gmm_fs_theta2) - GN(bhat_gmm_fs))/(0.001*bhat_gmm_fs(2, 1));
J_fs(1, 2) = dGNdtheta2_fs(1, 1);
J_fs(2, 2) = dGNdtheta2_fs(2, 1);
J_fs(3, 2) = dGNdtheta2_fs(3, 1);
clear bhat_gmm_fs_theta2 dGNdtheta2_fs;

% Obtaining the (1, 3), (2, 3), and (3, 3) components
bhat_gmm_fs_theta3 = bhat_gmm_fs;
bhat_gmm_fs_theta3(3, 1) = bhat_gmm_fs(3, 1)*1.001;
dGNdtheta3_fs = (GN(bhat_gmm_fs_theta3) - GN(bhat_gmm_fs))/(0.001*bhat_gmm_fs(3, 1));
J_fs(1, 3) = dGNdtheta3_fs(1, 1);
J_fs(2, 3) = dGNdtheta3_fs(2, 1);
J_fs(3, 3) = dGNdtheta3_fs(3, 1);
clear bhat_gmm_fs_theta3 dGNdtheta3_fs;
Gamma_fs = abs(J_fs);

% Obtaining GMM estimator var-cov matrix
varcov_matrix_gmm_fs = inv(Gamma_fs)*inv(A_bhat_gmm_fs)*inv(Gamma_fs');
% Obtaining numerical var and SE of (theta1, theta2, sigma)
var_bhat_gmm_fs = diag(varcov_matrix_gmm_fs);
clear varcov_matrix_gmm_fs J_fs Gamma_fs;
se_bhat_gmm_fs = sqrt(var_bhat_gmm_fs);
clear var_bhat_gmm_fs;

%==============================
% FIRST STAGE ESTIMATION END
% SECOND STAGE ESTIMATION START
%==============================

% Creating GMM objective function for second stage. This uses bhat_gmm_fs
% to recalculate the weight matrix A
objfunc_gmm = @(theta) GN(theta)'*A_bhat_gmm_fs*GN(theta);

% For second stage, we use bhat_gmm_fs as our initial value for the
% numerical optimisation.
theta0 = bhat_gmm_fs';
bhat_gmm = fminsearch(@(theta) objfunc_gmm(theta), theta0);
clear theta0;

% Calculating the A matrix using bhat_gmm found in the second stage.
% This will be used for creating the second stage varcov matrix.
A_bhat_gmm = A(bhat_gmm);

% Recall that GN is 3x1 vector. This means if we are taking the derivative
% of GN with respect to bhat_gmm', we are doing vector-by-vector
% derivatives. Our result will be a 3x3 matrix, or the Jacobian matrix.
% The Jacobian matrix is what we call Gamma/G at the bottom of this code.
% Creating Jacobian matrix
J = zeros(3, 3);

% Making bhat_gmm into a column vector
bhat_gmm = bhat_gmm';

% Obtaining the (1, 1), (2, 1), and (3, 1) components
bhat_gmm_theta1 = bhat_gmm;
bhat_gmm_theta1(1, 1) = bhat_gmm(1, 1)*1.001;
dGNdtheta1 = (GN(bhat_gmm_theta1) - GN(bhat_gmm))/(0.001*bhat_gmm(1, 1));
J(1, 1) = dGNdtheta1(1, 1);
J(2, 1) = dGNdtheta1(2, 1);
J(3, 1) = dGNdtheta1(3, 1);
clear bhat_gmm_theta1 dGNdtheta1;

% Obtaining the (1, 2), (2, 2), and (3, 2) components
bhat_gmm_theta2 = bhat_gmm;
bhat_gmm_theta2(2, 1) = bhat_gmm(2, 1)*1.001;
dGNdtheta2 = (GN(bhat_gmm_theta2) - GN(bhat_gmm))/(0.001*bhat_gmm(2, 1));
J(1, 2) = dGNdtheta2(1, 1);
J(2, 2) = dGNdtheta2(2, 1);
J(3, 2) = dGNdtheta2(3, 1);
clear bhat_gmm_theta2 dGNdtheta2;

% Obtaining the (1, 3), (2, 3), and (3, 3) components
bhat_gmm_theta3 = bhat_gmm;
bhat_gmm_theta3(3, 1) = bhat_gmm(3, 1)*1.001;
dGNdtheta3 = (GN(bhat_gmm_theta3) - GN(bhat_gmm))/(0.001*bhat_gmm(3, 1));
J(1, 3) = dGNdtheta3(1, 1);
J(2, 3) = dGNdtheta3(2, 1);
J(3, 3) = dGNdtheta3(3, 1);
clear bhat_gmm_theta3 dGNdtheta3;
Gamma = abs(J);

% Obtaining GMM estimator var-cov matrix
varcov_matrix_gmm = inv(Gamma)*inv(A_bhat_gmm)*inv(Gamma');
% Obtaining numerical var and SE of (theta1, theta2, sigma)
var_bhat_gmm = diag(varcov_matrix_gmm);
clear varcov_matrix_gmm J Gamma;
se_bhat_gmm = sqrt(var_bhat_gmm);
clear var_bhat_gmm;

%============================
% SECOND STAGE ESTIMATION END
%============================
%==========================================================================

%==========================================================================
%% Part of 3e that doesn't involve code

% We see that estimating the model without the (p_i)^2 moment results in
% the standard errors of all parameters to increase. This is expected
% because in general, GMM/MSM estimation does not use all of the moments of
% p_i's distribution (and if it does, it becomes extremely difficult to
% compute). In the context of this problem, this means the moment (p_i)^2
% was meant to capture some characteristics about the variance of p_i.
% Omitting it will therefore result in the GMM estimation to be less
% efficient.
%==========================================================================