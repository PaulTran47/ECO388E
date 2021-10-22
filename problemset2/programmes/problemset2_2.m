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

% ECO388E Problem Set 2, 2
% Paul Le Tran, plt377
% 20 October, 2021
%==========================================================================

%==========================================================================
%% Model info
% y_i = I(theta0 + theta1*x1_i + theta2*x2_i + e_i > 0);
% x2_i = theta3 + theta4*x1_i + theta5*z_i + eta_i
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
%% Part 2d: Joint normal discrete choice model (MLE): Estimation
% Assuming e_i and eta_i are jointly normal with SE(e_i) = 1 and SE(eta_i)
% = sigma_eta.
% First stage: x2_i = regressand; x1_i, z_i (plus constant) = regressors
% Second stage: y_i = regressand; x1_i, x2_i (plus constant) = regressors

%=====
% NOTE
%=====
% We are jointly estimating the two equations using MLE. To do this, we
% need to derive and calculate p(y_i, x2_i | x1_i, z_i; theta).
% For derivation, we first do this by making this joint probability the
% product of a marginal and conditional likelihood:

% p(y_i, x2_i | x1_i, z_i; theta) = p(x2_i | x1_i, z_i; theta)*p(y_i | x2_i, x1_i, z_i; theta)
%                                 = p(x2_i | x1_i, z_i; theta)*p(y_i | x2_i, x1_i, z_i, eta_i; theta).

% Consider the first probability on the RHS. We know x2_i and eta_i are
% invertible. For this reason, we can write:

% p(x2_i | x1_i, z_i; theta) = p(eta_i | x2_i, z_i; theta)*abs((dx2_i/deta_i))
%                            = normpdf(eta_i/sigma_eta)/sigma_eta.

% Consider the second probability on the RHS. We see that if:

% p(y_i == 1 | ...) = P(e_i + kappa_i < 0);
% p(y_i == 0 | ...) = 1 - P(e_i + kappa_i < 0),

% where kappa_i = theta0 + theta2*theta3 + (theta1 + theta2*theta4)*x1_i + theta2*theta5*z_i + theta2*eta_i.

% Given our assumption on the joint distribution of e_i and eta_i, we can
% derive the conditional distribution of e_i | eta_i as:

% e_i | eta_i ~ N((rho/(sigma_eta^2))*eta_i, 1 - (rho/sigma_eta)^2).

% Therefore, we can say:

% p(y_i | x2_i, x1_i, z_i, eta_i; theta) = (1 - y_i)*(1 - Phi(kappa_i)) + y_i*Phi(kappa_i),
% where Phi() is the above conditional CDF.
%=========
% END NOTE
%=========

%=====
% NOTE
%=====
% Mapping of parameter labels to model coefficients:
% theta1 = theta0;
% theta2 = theta1;
% theta3 = theta2;
% theta4 = theta3;
% theta5 = theta4;
% theta6 = theta5;
% theta7 = sigma_eta;
% theta8 = rho.
%=========
% END NOTE
%=========

% Creating inverse expression, eta_i
eta_i = @(theta) x2_i - theta(4) - theta(5)*x1_i - theta(6)*z_i;
% Creating p(eta_i | x2_i, z_i; theta)
p_eta_i = @(theta) normpdf(eta_i(theta)./abs(theta(7)));
% Creating p(x2_i | x1_i, z_i; theta)
p_x2_i = @(theta) p_eta_i(theta)./abs(theta(7));
% Creating kappa_i
kappa_i = @(theta) theta(1) + theta(3)*theta(4) + (theta(2) + theta(3)*theta(5))*x1_i + theta(3)*theta(6)*z_i + theta(3)*eta_i(theta);
% Creating mean of e_i | eta_i ~ N((rho/(sigma_eta^2))*eta_i, 1 - (rho/sigma_eta)^2)
mu_e_i = @(theta) (theta(8)/(abs(theta(7))^2))*eta_i(theta);
% Creating variance of e_i | eta_i ~ N((rho/(sigma_eta^2))*eta_i, 1 - (rho/sigma_eta)^2)
var_e_i = @(theta) 1 - (theta(8)/abs(theta(7)))^2;
% Creating p(y_i | x2_i, x1_i, z_i, eta_i; theta)
p_y_i = @(theta) (1 - y_i).*(1 - normcdf((kappa_i(theta) - mu_e_i(theta))./sqrt(var_e_i(theta)))) + y_i.*normcdf((kappa_i(theta) - mu_e_i(theta))./sqrt(var_e_i(theta)));
% Creating joint likelihood individual contribution of i, p(y_i, x2_i | x1_i, z_i; theta)
l_i = @(theta) p_x2_i(theta).*p_y_i(theta);
% Creating log-likelihood individual contribution of i
ll_i = @(theta) log(l_i(theta));
% Creating log-likelihood objective function
objfunc_ml = @(theta) sum(ll_i(theta))/N;

% Initial parameter values for numerical optimisation
theta0 = [-0.39, 0.01, 0.08, 9.1, 0.0, 0.36, 1.9, -0.11];

% Numerical optimisation to find (theta0, theta1, sigma)(fminsearch)
bhat_ml_joint = fminsearch(@(theta) -objfunc_ml(theta), theta0);
clear theta0;

% Calculating the estimated gradients of the log-likelihood contributions
% of i. We are estimating the gradient components numerically by a 
% "two-sided derivative" with steps of 0.001.

% Making parameter vector to column  vector
bhat_ml_joint = bhat_ml_joint';

% 1st gradient component for i
bhat_ml_joint_p_theta0 = bhat_ml_joint;
bhat_ml_joint_p_theta0(1, 1) = bhat_ml_joint(1, 1)*1.001;
bhat_ml_joint_m_theta0 = bhat_ml_joint;
bhat_ml_joint_m_theta0(1, 1) = bhat_ml_joint(1, 1)*0.999;
dll_i_jointdtheta0 = (ll_i(bhat_ml_joint_p_theta0) - ll_i(bhat_ml_joint_m_theta0))/(2*0.001*bhat_ml_joint(1, 1));
clear bhat_ml_joint_p_theta0 bhat_ml_joint_m_theta0;

% 2nd gradient component for i
bhat_ml_joint_p_theta1 = bhat_ml_joint;
bhat_ml_joint_p_theta1(2, 1) = bhat_ml_joint(2, 1)*1.001;
bhat_ml_joint_m_theta1 = bhat_ml_joint;
bhat_ml_joint_m_theta1(2, 1) = bhat_ml_joint(2, 1)*0.999;
dll_i_jointdtheta1 = (ll_i(bhat_ml_joint_p_theta1) - ll_i(bhat_ml_joint_m_theta1))/(2*0.001*bhat_ml_joint(2, 1));
clear bhat_ml_joint_p_theta1 bhat_ml_joint_m_theta1;

% 3nd gradient component for i
bhat_ml_joint_p_theta2 = bhat_ml_joint;
bhat_ml_joint_p_theta2(3, 1) = bhat_ml_joint(3, 1)*1.001;
bhat_ml_joint_m_theta2 = bhat_ml_joint;
bhat_ml_joint_m_theta2(3, 1) = bhat_ml_joint(3, 1)*0.999;
dll_i_jointdtheta2 = (ll_i(bhat_ml_joint_p_theta2) - ll_i(bhat_ml_joint_m_theta2))/(2*0.001*bhat_ml_joint(3, 1));
clear bhat_ml_joint_p_theta2 bhat_ml_joint_m_theta2;

% 4th gradient component for i
bhat_ml_joint_p_theta3 = bhat_ml_joint;
bhat_ml_joint_p_theta3(4, 1) = bhat_ml_joint(4, 1)*1.001;
bhat_ml_joint_m_theta3 = bhat_ml_joint;
bhat_ml_joint_m_theta3(4, 1) = bhat_ml_joint(4, 1)*0.999;
dll_i_jointdtheta3 = (ll_i(bhat_ml_joint_p_theta3) - ll_i(bhat_ml_joint_m_theta3))/(2*0.001*bhat_ml_joint(4, 1));
clear bhat_ml_joint_p_theta3 bhat_ml_joint_m_theta3;

% 5th gradient component for i
bhat_ml_joint_p_theta4 = bhat_ml_joint;
bhat_ml_joint_p_theta4(5, 1) = bhat_ml_joint(5, 1)*1.001;
bhat_ml_joint_m_theta4 = bhat_ml_joint;
bhat_ml_joint_m_theta4(5, 1) = bhat_ml_joint(5, 1)*0.999;
dll_i_jointdtheta4 = (ll_i(bhat_ml_joint_p_theta4) - ll_i(bhat_ml_joint_m_theta4))/(2*0.001*bhat_ml_joint(5, 1));
clear bhat_ml_joint_p_theta4 bhat_ml_joint_m_theta4;

% 6th gradient component for i
bhat_ml_joint_p_theta5 = bhat_ml_joint;
bhat_ml_joint_p_theta5(6, 1) = bhat_ml_joint(6, 1)*1.001;
bhat_ml_joint_m_theta5 = bhat_ml_joint;
bhat_ml_joint_m_theta5(6, 1) = bhat_ml_joint(6, 1)*0.999;
dll_i_jointdtheta5 = (ll_i(bhat_ml_joint_p_theta5) - ll_i(bhat_ml_joint_m_theta5))/(2*0.001*bhat_ml_joint(6, 1));
clear bhat_ml_joint_p_theta5 bhat_ml_joint_m_theta5;

% 7th gradient component for i
bhat_ml_joint_p_theta7 = bhat_ml_joint;
bhat_ml_joint_p_theta7(7, 1) = bhat_ml_joint(7, 1)*1.001;
bhat_ml_joint_m_theta7 = bhat_ml_joint;
bhat_ml_joint_m_theta7(7, 1) = bhat_ml_joint(7, 1)*0.999;
dll_i_jointdtheta7 = (ll_i(bhat_ml_joint_p_theta7) - ll_i(bhat_ml_joint_m_theta7))/(2*0.001*bhat_ml_joint(7, 1));
clear bhat_ml_joint_p_theta7 bhat_ml_joint_m_theta7;

% 8th gradient component for i
bhat_ml_joint_p_theta8 = bhat_ml_joint;
bhat_ml_joint_p_theta8(8, 1) = bhat_ml_joint(8, 1)*1.001;
bhat_ml_joint_m_theta8 = bhat_ml_joint;
bhat_ml_joint_m_theta8(8, 1) = bhat_ml_joint(8, 1)*0.999;
dll_i_jointdtheta8 = (ll_i(bhat_ml_joint_p_theta8) - ll_i(bhat_ml_joint_m_theta8))/(2*0.001*bhat_ml_joint(8, 1));
clear bhat_ml_joint_p_theta8 bhat_ml_joint_m_theta8;

% Creating numerical gradient for i
dll_i_jointdtheta = [dll_i_jointdtheta0 dll_i_jointdtheta1 dll_i_jointdtheta2 dll_i_jointdtheta3 dll_i_jointdtheta4 dll_i_jointdtheta5 dll_i_jointdtheta7 dll_i_jointdtheta8];
clear dll_i_jointdtheta0 dll_i_jointdtheta1 dll_i_jointdtheta2 dll_i_jointdtheta3 dll_i_jointdtheta4 dll_i_jointdtheta5 dll_i_jointdtheta7 dll_i_jointdtheta8;
dll_i_jointdtheta = dll_i_jointdtheta';

% Calculating numerical var-cov matrix
varcov_matrix_ml_joint = zeros(8, 8);
for i = 1:N
  varcov_matrix_ml_joint_i = dll_i_jointdtheta(:, i)*dll_i_jointdtheta(:, i)';
  varcov_matrix_ml_joint = varcov_matrix_ml_joint + varcov_matrix_ml_joint_i;
end
clear i varcov_matrix_ml_joint_i dll_i_jointdtheta;
varcov_matrix_ml_joint = inv(varcov_matrix_ml_joint);

% Obtaining numerical var and SE of (theta, sigma_eta, and rho)
var_bhat_ml_joint = diag(varcov_matrix_ml_joint);
clear varcov_matrix_ml_joint;
se_bhat_ml_joint = sqrt(var_bhat_ml_joint);
clear var_bhat_ml_joint;
%==========================================================================

%==========================================================================
%% Part 2e: Joint normal discrete choice model (MLE): Hypothesis-testing on rho
%=====
% NOTE
%=====
% We are testing if x2_i is correlated with e_i. In other words, we wish to
% see if the estimated value of rho is statistically significant. We will
% do this with a simple two-tailed t-test and alpha = 0.05. Our null
% hypothesis is that rho = 0. Finally, df = N - 3 because we have three
% variables: x1_i, x2_i, and z_i.
%=========
% END NOTE
%=========
tstat = bhat_ml_joint(8, 1)/(se_bhat_ml_joint(8,1)/sqrt(N));
p_value = tcdf(tstat, N - 3);

%=======
% ANSWER
%=======
% Given our t-statistic (with df = 752) and associated p-value, we reject
% the null hypothesis at an alpha level of 0.05. As a result, we can say
% it's possible that x2_i is correlated with e_i, and thus, x2_i is
% endogenous.
%===========
% END ANSWER
%===========
%==========================================================================
%% Questions asked in problemset 2_2 that don't involve code
%=========
% Part 2a:
%=========
% Recall that x2_i is women's education level (in years). One economic
% reason how x2_i could be correlated with e_i is that e_i could be
% capturing the wage a woman wants from working. This is because women who
% want a higher wage may choose more education. Therefore, if women make
% this decision simultaneously regarding the two variables, you would
% expect the correlation between x2_i and e_i to put upward bias on
% coefficient theta2.

%=========
% Part 2b:
%=========
% NO QUESTION TO ANSWER. THIS PART IS SETUP FOR PART 2c.

%=========
% Part 2c:
%=========
% From part 2b, we assume that a woman's educational level depends on x1_i,
% z_i, and e_i. From the system of two equations, we assume that
% unobservables e_i and eta_i are jointly normal and are independent of
% x1_i and z_i. Observe that
% cov(e_i, eta_i) = rho = E[eta_i*e_i] - E[eta_i]*E[e_i] = E[eta_i*e_i],
% by the jointly normal assumption.

% Rho relates to the possible endogeneity of x2_i through expectations.
% Specifically, with the given system of two equations, we know that
% E[x2_i*e_i] = E[eta_i*e_i] = rho. Therefore, we can see that if x2_i is
% exogenous, we must have rho = 0. But since x2_i is assumed to be
% endogenous, rho != 0. Additionally, I expect rho > 0. This is because the
% educational level of one's parents (z_i) should be positively correlated
% with x2_i. This is because a parent often serves as a role model.

% Finally, the assumption that z_i does not directly determine y_i
% (i.e., z_i affects y_i only through x2_i) is exclusion restriction. This
% is one of the three principles that makes for a "good"/"robust"
% instrument in IV estimation. Without the exclusion restriction, theta2
% would not be identified non-parametrically.

%=========
% Part 2f:
%=========
% Observe that tau_i is an unobservable that affects the effect of x1_i
% (women's age) on x2_i randomly. Economically, this means tau_i is
% capturing heterogeneity in the response of education from changes in age.
% In summary, a woman's age could have a larger impact on her educational
% level compared to other women.

%=========
% Part 2g:
%=========
% Note that we are choosing to partially integrate out tau_i.

% Similar to the model in part 2d, we can write the desired
% joint likelihood as:

% p(y_i, x2_i | x1_i, z_i; theta) = p(x2_i | x1_i, z_i; theta)*p(y_i | x1_i, z_i; theta)
%                                 = p(x2_i | x1_i, z_i; theta)*p(y_i | x1_i, x2_i, z_i, eta_i; theta).

% The last equality holds because like before, the first stage equation is
% still invertible. In other words, this holds because
% p(x2_i | x1_i, z_i; theta) = p(x2_i | x1_i, z_i, eta_i; theta).
% This means even with the additional unobservable tau_i, we do not need to
% integrate out to find p(x2_i | x1_i, z_i; theta).

% For p(x2_i | x1_i, z_i; theta), x2_i is still invertible in eta_i. As a
% result, we can write the following expression:
% p(x2_i | x1_i, z_i; theta) = p(x2_i | x1_i, z_i, eta_i; theta)
% However, we see that the RHS can only be found if we partially integrate
% out an unobservable. As mentioned in the beginning, we choose tau_i:

% p(x2_i | x1_i, z_i, eta_i; theta) = \int p(x2_i | x1_i, z_i, eta_i, tau_i; theta)*p(tau_i)dtau_i.

% We are able to use the unconditional p(tau_i) because of assumption that 
% tau_i is independent of x1_i, x2_i, z_i, eta_i, and e_i. Finally, because
% x2_i is invertible in eta_i, we apply the ToRV:

% p(x2_i | x1_i, z_i, eta_i; theta) = \int p(eta_i | x1_i, z_i, x2_i, tau_i; theta)*abs(det(deta_i/dy_i))*p(tau_i)dtau_i.
%                                   = inv(S)*\sum_{j}^{S} normpdf(eta_i/sigma_eta)/sigma_eta,

% where
% eta_i = x2_i - theta3 - (theta4 + sigma_tau*tau_{i, S})x1_i - theta5*z_i.
% S represents the number of simulated draws we take of tau_i and
% j = 1, ..., S is simulated draw itself.
% We can finally write the probabilty as:

% p(x2_i | x1_i, z_i; theta) = inv(S)*\sum_{j}^{S} normpdf(x2_i - theta3 - (theta4 + sigma_tau*tau_{i, j})x1_i - theta5*z_i)/sigma_eta

% To find an expression for for p(y_i | x1_i, x2_i, z_i, eta_i; theta), we
% again partially integrate out tau_i:

% p(y_i | x1_i, x2_i, z_i, eta_i; theta) = \int p(y_i | x1_i, x2_i, z_i, eta_i, tau_i; theta)*p(tau_i)dtau_i.

% We are able to use the unconditional p(tau_i) because of assumption that 
% tau_i is independent of x1_i, x2_i, z_i, eta_i, and e_i.
% Simulating this integral gives us:

% p(y_i | x1_i, x2_i, z_i, eta_i; theta) = inv(S)*\sum_{j}^{S} p(y_i | x1_i, x2_i, z_i, tau_{i, j}, eta_i; theta),

% We must now modify our kappa_i equation that is in the second stage
% equation due to the structure of x2_i changing. Essentially, this just
% means kappa_i now includes tau_i:

% kappa_i = theta0 + theta2*theta3 + (theta1 + theta2*(theta4 + sigma_tau*tau_{i, j}))*x1_i + theta2*theta5*z_i + theta2*eta_i.

% Therefore, we can write the our simulated (log-)likelihood objective
% function:

% L = \sum_{i}^{N} ln(p(x2_i | x1_i, z_i; theta)) + ln(p(y_i | x1_i, x2_i, z_i, eta_i; theta)),

% where
% p(y_i | x1_i, x2_i, z_i, eta_i; theta) = inv(S)*\sum_{j}^{S} y_i*normcdf((kappa_i + (rho/sigma_eta)*eta_i)/(1 - rho^2)) + (1 - y_i)*normcdf((kappa_i + (rho/sigma_eta)*eta_i)/(1 - rho^2)),
% where
% eta_i = x2_i - theta3 - (theta4 + sigma_tau*tau_{i, j})x1_i - theta5*z_i,
% kappa_i = theta0 + theta2*theta3 + (theta1 + theta2*(theta4 + sigma_tau*tau_{i, j}))*x1_i + theta2*theta5*z_i + theta2*eta_i,
% S represents the number of simulated draws we take of tau_i,
% j = 1, ..., S is simulated draw itself.

%=========
% Part 2h:
%=========
% Observe that we now have 9 parameters in our model: theta0 - theta5,
% rho, sigma_eta, and sigma_tau. We aim to choose 9 moments for MSM
% estimation. Furthermore, observe that our model is basically an
% application of the general case of endogeneity in non-invertible models
% section of class. Because of this, it is simplest to form generic moments
% using the "reduced form" of our model:

% Let \widetilde{y_i} = (y_i, x2_i), \widetilde{e_i} = (e_i, eta_i, tau_i).

% Then we can write the following "generic" moments:
% i. \widetilde{y_i} - E[\widetilde{y_i} | x1_i, z_i; theta];
% ii. (\widetilde{y_i} - E[\widetilde{y_i} | x1_i, z_i; theta])*x1_i;
% iii. (\widetilde{y_i} - E[\widetilde{y_i} | x1_i, z_i; theta])*z_i;
% iv. \widetilde{y_i}^2 - E[\widetilde{y_i}^2 | x1_i, z_i; theta];

% Translating the "generic" moments of the "reduced form" to the original
% model yields the following eight moments:
% 1. y_i - E[y_i | x1_i, z_i; theta];
% 2. (y_i - E[y_i | x1_i, z_i; theta])*x1_i;
% 3. (y_i - E[y_i | x1_i, z_i; theta])*z_i;
% 4. x2_i - E[x2_i | x1_i, z_i; theta];
% 5. (x2_i - E[x2_i | x1_i, z_i; theta])*x1_i;
% 6. (x2_i - E[x2_i | x1_i, z_i; theta])*z_i;
% 7. y_i^2 - E[y_i^2 | x1_i, z_i; theta];
% 8. x2_i^2 - E[x2_i^2 | x1_i, z_i; theta];

% Moments seven and eight (iv) will help in identifying parameters rho and
% sigma_eta. This is because the squared terms capture more of the
% endogeneity effect. Recall our independence assumption about tau_i. From
% it, we can state that

% E[tau_i*x1_i] = E[tau_i*z_i] = 0;

% Therefore, using the idea that squared terms capture more information
% about variance related parameters, the moments nine and 10 should help
% identify sigma_tau:
% 9. E[tau_i*x1_i^2] = 0;
% 10. E[tau_i*z_i^2] = 0.

%=========
% Part 2i:
%=========
% The SML estimator has several advantages over the MSM estimator. For one,
% the SML estimator is always efficient in that all the necessary moments
% are included in the estimation. In contrast, the MSM estimator is only
% efficient if we choose the optimal moments function. Otherwise, it will
% not be so. Additionally, the SML estimator in part 2g has the advantange
% in that we only need to simulate draws for tau_i.

% The MSM estimator has the big advantage in that it will remain
% consistent holding the number of simulated draws we take of the error
% terms integrated out constant. In contrast, the SML estimator requires
% the set of simulated draws to approach infinity in order to be
% consistent.

%=========
% Part 2j:
%=========
% Soon TM.

%=========
% Part 2k:
%=========
% It's possible that assuming tau_i is independent of eta_i is not
% reasonable. Recall from part 2f that tau_i captures heterogeneity in the
% response of a woman's education to changes in her age. Given the
% endogenous nature of our model, if education and labour force
% participation are determined at the same time by a woman, it's possible
% that eta_i captures a factor in that simultaneous decision. It's
% additionally possible that the "shock" affecting a woman's choice of
% education due to her age has influence on this simultaneous decision on
% education and labour force participation as well. As a result, assuming
% tau_i is independent of eta_i could be unreasonable.
%==========================================================================