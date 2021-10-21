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

% ECO388E Problem Set 2, 1
% Paul Le Tran, plt377
% 20 October, 2021
%==========================================================================

%==========================================================================
%% Model info
% y_i = I(theta1 + theta2*x1_i + theta3*x2_i + e_i > 0)
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
data1 = importdata('ps2.dat');

% Declaring variables
y_i = data1(:, 1);
x1_i = data1(:, 2);
x2_i = data1(:, 3);
z_i = data1(:, 4);

% Getting sample size
% Variable is global so it can be accessed in the separate workspaces of
% called functions
global N;
N = length(y_i);

clear data1;
cd(home_dir);
%==========================================================================

%==========================================================================
%% Part 1a: Binary probit model (MLE): Estimation
% Assuming e_i ~ N(0, 1), iid
% y_i = regressand; x1_i, x2_i (plus constant) = regressors
% Recall that the likelihood contribution of i in a binary probit model is
% p(y_i = 1 | x1_i, x2_i; theta1, theta2, theta3) = normcdf(theta1 + theta2*x1_i + theta3*x2_i)

%=============
% MLE manually
%=============
% Creating p(y_i = 1 | x1_i, x2_i; theta1, theta2, theta3)
l_i_bp = @(theta) normcdf(theta(1) + theta(2)*x1_i + theta(3)*x2_i);
% Creating log-likelihood contribution of i
ll_i_bp = @(theta) log(l_i_bp(theta));
% Creating log-likelihood objective function
objfunc_ml_bp = @(theta) sum(y_i.*ll_i_bp(theta) + (1 - y_i).*log(1 - l_i_bp(theta)))/N;

% Initial parameter values for numerical optimisation
theta0 = [0, 0, 0];

% Numerical optimisation to find (theta1, theta2, sigma)(fminsearch)
bhat_ml_bl = fminsearch(@(theta) -objfunc_ml_bp(theta), theta0);
clear theta0;

% Calculating the estimated gradients of the log-likelihood contributions
% of i. We are estimating the gradient components numerically by a 
% "two-sided derivative" with steps of 0.001.

% Making parameter vector to column  vector
bhat_ml_bl = bhat_ml_bl';

% 1st gradient component for i
bhat_ml_bl_p_theta1 = bhat_ml_bl;
bhat_ml_bl_p_theta1(1, 1) = bhat_ml_bl(1, 1)*1.001;
bhat_ml_bl_m_theta1 = bhat_ml_bl;
bhat_ml_bl_m_theta1(1, 1) = bhat_ml_bl(1, 1)*0.999;
dll_i_bldtheta1 = (ll_i_bp(bhat_ml_bl_p_theta1) - ll_i_bp(bhat_ml_bl_m_theta1))/(2*0.001*bhat_ml_bl(1, 1));
clear bhat_ml_bp_p_theta1 bhat_ml_bp_m_theta1;

% 2nd gradient component for i
bhat_ml_bl_p_theta2 = bhat_ml_bl;
bhat_ml_bl_p_theta2(2, 1) = bhat_ml_bl(2, 1)*1.001;
bhat_ml_bl_m_theta2 = bhat_ml_bl;
bhat_ml_bl_m_theta2(2, 1) = bhat_ml_bl(2, 1)*0.999;
dll_i_bldtheta2 = (ll_i_bp(bhat_ml_bl_p_theta2) - ll_i_bp(bhat_ml_bl_m_theta2))/(2*0.001*bhat_ml_bl(2, 1));
clear bhat_ml_bp_p_theta2 bhat_ml_bp_m_theta2;

% 3nd gradient component for i
bhat_ml_bl_p_theta3 = bhat_ml_bl;
bhat_ml_bl_p_theta3(3, 1) = bhat_ml_bl(3, 1)*1.001;
bhat_ml_bl_m_theta3 = bhat_ml_bl;
bhat_ml_bl_m_theta3(3, 1) = bhat_ml_bl(3, 1)*0.999;
dll_i_bldtheta3 = (ll_i_bp(bhat_ml_bl_p_theta3) - ll_i_bp(bhat_ml_bl_m_theta3))/(2*0.001*bhat_ml_bl(3, 1));
clear bhat_ml_bp_p_theta3 bhat_ml_bp_m_theta3;

% Creating numerical gradient for i
dll_i_bldtheta = [dll_i_bldtheta1 dll_i_bldtheta2 dll_i_bldtheta3];
clear dll_i_bpdtheta1 dll_i_bpdtheta2 dll_i_bpdtheta3;
dll_i_bldtheta = dll_i_bldtheta';

% Calculating numerical var-cov matrix
varcov_matrix_ml_bl = zeros(3, 3);
for i = 1:N
  varcov_matrix_ml_bl_i = dll_i_bldtheta(:, i)*dll_i_bldtheta(:, i)';
  varcov_matrix_ml_bl = varcov_matrix_ml_bl + varcov_matrix_ml_bl_i;
end
clear i varcov_matrix_ml_bp_i dll_i_bpdtheta;
varcov_matrix_ml_bl = inv(varcov_matrix_ml_bl);

% Obtaining numerical var and SE of (theta1, theta2, theta3)
var_bhat_ml_bl = diag(varcov_matrix_ml_bl);
clear varcov_matrix_ml_bp;
se_bhat_ml_bl = sqrt(var_bhat_ml_bl);
clear var_bhat_ml_bp;

%===============
% MLE via mnrfit
%===============
% We are not adding columns of ones in regressor matrix because mnrfit
% automatically includes a constant in the regression. We also are
% converting y_i into categorical class since mnrfit doesn't like doubles.
%=====
% NOTE
%=====
% For some reason, the signs are correct only when you use -y_i instead of
% y_i.
%=========
% END NOTE
%=========
[bhat_ml_bp_mnrfit, dev, stats] = mnrfit([x1_i x2_i], categorical(-y_i), 'model', 'ordinal', 'link', 'probit');
% Saving SEs
se_bhat_ml_bp_mnrfit = stats.se;
clear dev;
%=====
% NOTE
%=====
% The manual and mnrfit SEs are different, but by a small magnitude ~(0.03,
% 0.001, and 0.001). For this reason, I report both (since I'm confident
% the manual code is correct).
%=========
% END NOTE
%=========
%==========================================================================

%==========================================================================
%% Part 1b: Binary probit model (MLE): Marginal effects
% To have a better interpretation of bhat_ml, we are calculating the
% marginal effects using the average of an average woman.

% Saving averages of regressors into a vector
x_bar = [1 mean(x1_i) mean(x2_i)];

% Calculating average p(y_i = 1 | x1_i, x2_i; theta1, theta2, theta3) using
% bhat_ml
l_bar = normpdf(x_bar*bhat_ml_bl)*bhat_ml_bl;
clear x_bar;
%=======
% ANSWER
%=======
% We see that for a woman of mean age and education, the effect of an
% additional year of education on the probability of her working is roughly
% 4.1%.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 1c: Binary logit model (MLE)
% Assuming e_i ~ logistic (i.e., the difference between two iid extreme
% value deviates)
% p(y_i = 1 | x1_i, x2_i; theta1, theta2, theta3) = exp(theta1 +theta2*x1_i + theta3*x2_i)/(1 + exp(theta1 +theta2*x1_i + theta3*x2_i))

%=============
% MLE manually
%=============
% Creating p(y_i = 1 | x1_i, x2_i; theta1, theta2, theta3)
l_i_bl = @(theta) exp(theta(1) + theta(2)*x1_i + theta(3)*x2_i)./(1 + exp(theta(1) + theta(2)*x1_i + theta(3)*x2_i));
% Creating log-likelihood contribution of i
ll_i_bl = @(theta) log(l_i_bl(theta));
% Creating log-likelihood objective function
objfunc_ml_bl = @(theta) inv(N)*sum(y_i.*ll_i_bl(theta) + (1 - y_i).*log(1 - l_i_bl(theta)));
% Initial parameter values for numerical optimisation
theta0 = [0, 0, 0];

% Numerical optimisation to find (theta1, theta2, sigma)(fminsearch)
bhat_ml_bl = fminsearch(@(theta) -objfunc_ml_bl(theta), theta0);
clear theta0;

% Calculating the estimated gradients of the log-likelihood contributions
% of i. We are estimating the gradient components numerically by a 
% "two-sided derivative" with steps of 0.001.

% Making parameter vector to column  vector
bhat_ml_bl = bhat_ml_bl';

% 1st gradient component for i
bhat_ml_bl_p_theta1 = bhat_ml_bl;
bhat_ml_bl_p_theta1(1, 1) = bhat_ml_bl(1, 1)*1.001;
bhat_ml_bl_m_theta1 = bhat_ml_bl;
bhat_ml_bl_m_theta1(1, 1) = bhat_ml_bl(1, 1)*0.999;
dll_i_bldtheta1 = (ll_i_bl(bhat_ml_bl_p_theta1) - ll_i_bl(bhat_ml_bl_m_theta1))/(2*0.001*bhat_ml_bl(1, 1));
clear bhat_ml_bl_p_theta1 bhat_ml_bl_m_theta1;

% 2nd gradient component for i
bhat_ml_bl_p_theta2 = bhat_ml_bl;
bhat_ml_bl_p_theta2(2, 1) = bhat_ml_bl(2, 1)*1.001;
bhat_ml_bl_m_theta2 = bhat_ml_bl;
bhat_ml_bl_m_theta2(2, 1) = bhat_ml_bl(2, 1)*0.999;
dll_i_bldtheta2 = (ll_i_bl(bhat_ml_bl_p_theta2) - ll_i_bl(bhat_ml_bl_m_theta2))/(2*0.001*bhat_ml_bl(2, 1));
clear bhat_ml_bl_p_theta2 bhat_ml_bl_m_theta2;

% 3nd gradient component for i
bhat_ml_bl_p_theta3 = bhat_ml_bl;
bhat_ml_bl_p_theta3(3, 1) = bhat_ml_bl(3, 1)*1.001;
bhat_ml_bl_m_theta3 = bhat_ml_bl;
bhat_ml_bl_m_theta3(3, 1) = bhat_ml_bl(3, 1)*0.999;
dll_i_bldtheta3 = (ll_i_bl(bhat_ml_bl_p_theta3) - ll_i_bl(bhat_ml_bl_m_theta3))/(2*0.001*bhat_ml_bl(3, 1));
clear bhat_ml_bl_p_theta3 bhat_ml_bl_m_theta3;

% Creating numerical gradient for i
dll_i_bldtheta = [dll_i_bldtheta1 dll_i_bldtheta2 dll_i_bldtheta3];
clear dll_i_bldtheta1 dll_i_bldtheta2 dll_i_bldtheta3;
dll_i_bldtheta = dll_i_bldtheta';

% Calculating numerical var-cov matrix
varcov_matrix_ml_bl = zeros(3, 3);
for i = 1:N
  varcov_matrix_ml_bl_i = dll_i_bldtheta(:, i)*dll_i_bldtheta(:, i)';
  varcov_matrix_ml_bl = varcov_matrix_ml_bl + varcov_matrix_ml_bl_i;
end
clear i varcov_matrix_ml_bl_i dll_i_bldtheta;
varcov_matrix_ml_bl = inv(varcov_matrix_ml_bl);

% Obtaining numerical var and SE of (theta1, theta2, theta3)
var_bhat_ml_bl = diag(varcov_matrix_ml_bl);
clear varcov_matrix_ml_bl;
se_bhat_ml_bl = sqrt(var_bhat_ml_bl);
clear var_bhat_ml_bl;

%===============
% MLE via mnrfit
%===============
% We are not adding columns of ones in regressor matrix because mnrfit
% automatically includes a constant in the regression. We also are
% converting y_i into categorical class since mnrfit doesn't like doubles.
%=====
% NOTE
%=====
% For some reason, the signs are correct only when you use -y_i instead of
% y_i.
%=========
% END NOTE
%=========
[bhat_ml_bl_mnrfit, dev, stats] = mnrfit([x1_i x2_i], categorical(-y_i), 'model', 'ordinal', 'link', 'logit');
% Saving SEs
se_bhat_ml_bl_mnrfit = stats.se;

%=======
% ANSWER
%=======
% Observe that the marginal effect based on the average of the data for
% both models has the same analytical structure: The marginal probability
% effect of regressor xj = normpdf(X'theta)*thetaj. Recall that the tails
% of the logit distribution are "fatter" than those of the probit
% distribution. As a result, the marginal effect of regressor xj under
% logit will be greater than that under probit because normpdf() under
% logit is larger in magnitude. This ultimately results in the estimated
% coefficients under logit to be larger in magnitude than those estimated
% under probit.

% However, the overall curvature of both distributions are similar. This is
% because the structure of the model (i.e., the relation between the
% variables and and the error term) has not changed. As a result, the ratio
% amongst coefficients within probit and logit don't change much at all.
%===========
% END ANSWER
%===========
%==========================================================================