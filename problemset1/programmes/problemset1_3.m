%==========================================================================
%% 2 percentage signs represent sections of code;
% 1 percentage sign represents comments for code or commented out code;

% ECO388E Problem Set 1, 3
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
data_drawsml = importdata('drawsml.dat');
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
%% Part 3b: ML
% Assuming e1_i, e2_i ~ N(0, sigma^2), iid; sigma > 0
% ==> sigma*e1_i, sigma*e2_i ~ N(0, 1), iid
% Therefore, model:
% p_i = (100 + exp(theta1 + theta2*x1_i + sigma*e1_i) + exp(theta1 + theta2*x2_i + sigma*e2_i))/3
% Observe that this model is both non-linear and non-invertible because
% we have more unobservables than number of regressands. Therefore, we need
% to integrate out one unobservable in order to do ML estimation.
%=======================================
% WE ARE CHOOSING TO INTEGRATE OUT e1_i:
%=======================================
% Doing so shows that there is an invertible relationship between p_i and
% e2_i, we can do ToRV:
% e2_i = (ln(3*p_i - 100 - exp(theta1 + theta2*x1_i + sigma*e1_i)) - theta1 - theta2*x2_i)/sigma
% d(f^(-1))/d(p_i) = 3/(3*p_i - 100 - exp(theta1 + theta2*x1_i + sigma*e1_i))
% Therefore, p(p_i|x1_i, x2_i, e1_i; theta1, theta2, sigma) =
% normpdf(e2_i)*abs(3/(3*p_i - 100 - exp(theta1 + theta2*x1_i + sigma*e1_i)))
% Therefore, using numerical approximation of integrals, we get the
% likelihood p(p_i|x1_i, x2_i, e1_i, e2_i; theta1, theta2, sigma) = 
% inv(S)*sum_over_S(p(p_i|x1_i, x2_i, e1_S_i; theta1, theta2, sigma)),
% where S represents the simulated draw.
% What this means is we need to calculate p(p_i|x1_i, x2_i, e1_S_i; theta1, theta2, sigma)
% for each simulated draw of e1_i values first.

% Creating matrix that houses all the simulated draws for e1_i (20 draws,
% so the matrix is 50x20 due to N = 50).
% Matrix is global so it can be accessed in the separate workspaces of
% called functions
global e1_i_matrix;
e1_i_matrix = data_drawsml;
% Creating variable storing number of simulated draws
% Variable is global so it can be accessed in the separate workspaces of
% called functions
global S;
S = min(size(e1_i_matrix));
clear data_drawsml;

% Creating p(p_i|x1_i, x2_i, e1_i, e2_i; theta1, theta2, sigma) = likelihood contribution of i
% More info about the l_io_e1_i_i_sum function can be found in
% l_io_e1_i_i_sum.m
l_i = @(theta) inv(S)*l_io_e1_i_i_sum(theta);
% Creating log-likelihood contribution of i
ll_i = @(theta) log(l_i(theta));
% Creating log-likelihood function
ll = @(theta) inv(N)*sum(ll_i(theta));

% Because we have a non-linear model, we need to test out different initial
% values to see what the global max is. For simplicity, we will restrict
% the domains of the parameters to be ([0, 1.5], [0, 1.5] [0, 1.5]). We
% will search for the highest local max via a loop. We will be skipping by
% 0.3 for every parameter to save time.
% Creating counter variable
counter = 0;
% Starting value of the initial parameter values
theta0_start = [0, 0, 0];
% Creating a matrix that will save every iteration's local max.
bhat_ml_matrix = zeros(216, 3);
for i = theta0_start(1, 1):0.3:1.5 % theta1 values
  for j = theta0_start(1, 2):0.3:1.5 %theta2 values
    for k = theta0_start(1, 3):0.3:1.5 %theta3/sigma values
      theta0 = [i, j, k]; % Actual initial parameter values for numerical optimisation
      % Numerical optimisation to find (theta1, theta2, sigma)(fminsearch)
      disp(counter);
      disp(theta0);
      bhat_ml = fminsearch(@(theta) -ll(theta), theta0);
      disp(bhat_ml);
      counter = counter + 1;
      bhat_ml_matrix(counter, :) = bhat_ml;
    end
  end
end
clear counter theta0_start i j k theta0 bhat_ml;

% Creating a vector that stores the log-likelihood when given each bhat_ml,
% and appending it as a column vector to the right of bhat_ml_matrix
ll_vector = zeros(length(bhat_ml_matrix), 1);
for i = 1:length(bhat_ml_matrix)
  ll_vector(i, 1) = ll(bhat_ml_matrix(i, :));
end
ll_bhat_ml_matrix = [ll_vector bhat_ml_matrix];

% We will be removing any instances of the log-likelihood being -Inf
% (which corresponds with integer optimal parameter values; these don't
% make sense when it comes to numerical optimisation).
condition1 = ll_bhat_ml_matrix(:, 1) == -Inf;
ll_bhat_ml_matrix(condition1, :) = [];
clear ll_vector bhat_ml_matrix condition1;

% Finding the bhat_ml values which give the maximum log-likelihood.
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
clear i dll_idtheta varcov_matrix_ml_i;
varcov_matrix_ml = inv(varcov_matrix_ml);
% Obtaining numerical var and SE of (theta1, theta2, sigma)
var_bhat_ml = diag(varcov_matrix_ml);
clear varcov_matrix_ml;
se_bhat_ml = sqrt(var_bhat_ml);
se_bhat_ml = se_bhat_ml(1:2, 1);
clear var_bhat_ml;
%==========================================================================

%==========================================================================
%% Part 3d: GMM (Technically MSM)
% Creating matrix that houses all the simulated draws for e1_i and e2_i
% Both have 20 draws, so the matrix is 50x40 due to N = 50. The first 20
% columns are draws for e1_i. The second 20 draws are draws for e2_i.
% Matrix is global so it can be accessed in the separate workspaces of
% called functions
global e_i_matrix;
e_i_matrix = data_drawsgmm;
clear data_drawsgmm;

% Creating the left component of the first moment in the individual moment
% function g_i. This is the difference between p_i and the simulated inner
% expectations of p_i (50x1). More info about the inner_expectations_p_sum
% function can be found in inner_expectations_p_sum.m
g_1st_m_component_i = @(theta) regressand - inv(S)*inner_expectations_p_sum(theta);

% Creating the individual first moment function g_1st_m_i (50x3)
g_1st_m_i = @(theta) [g_1st_m_component_i(theta) g_1st_m_component_i(theta).*x1_i g_1st_m_component_i(theta).*x2_i];

% Creating the left component of the second moment in the individual moment
% function g_i. This is the difference between (p_i)^2 and the simulated
% inner expectations of (p_i)^2 (50x1). More info about the
% inner_expectations_p2_sum function can be found in
% inner_expectations_p2_sum.m
g_2nd_m_component_i = @(theta) regressand.^2 - inv(S)*inner_expectations_p2_sum(theta);

% Creating the individual second moment function g_2nd_m_i (50x1)
g_2nd_m_i = @(theta) g_2nd_m_component_i(theta);

% Creating the individual moment function g_i (50x4)
g_i = @(theta) [g_1st_m_i(theta) g_2nd_m_i(theta)];

% Creating mean aggregate moment function GN_i (4x1)
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

% Recall that GN is 4x1 vector. This means if we are taking the derivative
% of GN with respect to bhat_gmm', we are doing vector-by-vector
% derivatives. Our result will be a 4x3 matrix, or the Jacobian matrix.
% The Jacobian matrix is what we call Gamma/G at the bottom of this code.
% Creating Jacobian matrix
J_fs = zeros(4, 3);

% Making bhat_gmm_fs into a column vector
bhat_gmm_fs = bhat_gmm_fs';

% Obtaining the (1, 1), (2, 1), (3, 1), and (4, 1) components
bhat_gmm_fs_theta1 = bhat_gmm_fs;
bhat_gmm_fs_theta1(1, 1) = bhat_gmm_fs(1, 1)*1.001;
dGNdtheta1_fs = (GN(bhat_gmm_fs_theta1) - GN(bhat_gmm_fs))/(0.001*bhat_gmm_fs(1, 1));
J_fs(1, 1) = dGNdtheta1_fs(1, 1);
J_fs(2, 1) = dGNdtheta1_fs(2, 1);
J_fs(3, 1) = dGNdtheta1_fs(3, 1);
J_fs(4, 1) = dGNdtheta1_fs(4, 1);
clear bhat_gmm_fs_theta1 dGNdtheta1_fs;

% Obtaining the (1, 2), (2, 2), (3, 2), and (4, 2) components
bhat_gmm_fs_theta2 = bhat_gmm_fs;
bhat_gmm_fs_theta2(2, 1) = bhat_gmm_fs(2, 1)*1.001;
dGNdtheta2_fs = (GN(bhat_gmm_fs_theta2) - GN(bhat_gmm_fs))/(0.001*bhat_gmm_fs(2, 1));
J_fs(1, 2) = dGNdtheta2_fs(1, 1);
J_fs(2, 2) = dGNdtheta2_fs(2, 1);
J_fs(3, 2) = dGNdtheta2_fs(3, 1);
J_fs(4, 2) = dGNdtheta2_fs(4, 1);
clear bhat_gmm_fs_theta2 dGNdtheta2_fs;

% Obtaining the (1, 3), (2, 3), (3, 3), and (4, 3) components
bhat_gmm_fs_theta3 = bhat_gmm_fs;
bhat_gmm_fs_theta3(3, 1) = bhat_gmm_fs(3, 1)*1.001;
dGNdtheta3_fs = (GN(bhat_gmm_fs_theta3) - GN(bhat_gmm_fs))/(0.001*bhat_gmm_fs(3, 1));
J_fs(1, 3) = dGNdtheta3_fs(1, 1);
J_fs(2, 3) = dGNdtheta3_fs(2, 1);
J_fs(3, 3) = dGNdtheta3_fs(3, 1);
J_fs(4, 3) = dGNdtheta3_fs(4, 1);
clear bhat_gmm_fs_theta3 dGNdtheta3_fs;
Gamma_fs = abs(J_fs);

% Obtaining GMM estimator var-cov matrix
varcov_matrix_gmm_fs = inv(Gamma_fs'*A_bhat_gmm_fs*Gamma_fs)*(Gamma_fs'*A_bhat_gmm_fs*inv(A_bhat_gmm_fs)*A_bhat_gmm_fs*Gamma_fs)*inv(Gamma_fs'*A_bhat_gmm_fs*Gamma_fs);
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

% Recall that GN is 4x1 vector. This means if we are taking the derivative
% of GN with respect to bhat_gmm', we are doing vector-by-vector
% derivatives. Our result will be a 4x3 matrix, or the Jacobian matrix.
% The Jacobian matrix is what we call Gamma/G at the bottom of this code.
% Creating Jacobian matrix
J = zeros(4, 3);

% Making bhat_gmm into a column vector
bhat_gmm = bhat_gmm';

% Obtaining the (1, 1), (2, 1), (3, 1), and (4, 1) components
bhat_gmm_theta1 = bhat_gmm;
bhat_gmm_theta1(1, 1) = bhat_gmm(1, 1)*1.001;
dGNdtheta1 = (GN(bhat_gmm_theta1) - GN(bhat_gmm))/(0.001*bhat_gmm(1, 1));
J(1, 1) = dGNdtheta1(1, 1);
J(2, 1) = dGNdtheta1(2, 1);
J(3, 1) = dGNdtheta1(3, 1);
J(4, 1) = dGNdtheta1(4, 1);
clear bhat_gmm_theta1 dGNdtheta1;

% Obtaining the (1, 2), (2, 2), (3, 2), and (4, 2) components
bhat_gmm_theta2 = bhat_gmm;
bhat_gmm_theta2(2, 1) = bhat_gmm(2, 1)*1.001;
dGNdtheta2 = (GN(bhat_gmm_theta2) - GN(bhat_gmm))/(0.001*bhat_gmm(2, 1));
J(1, 2) = dGNdtheta2(1, 1);
J(2, 2) = dGNdtheta2(2, 1);
J(3, 2) = dGNdtheta2(3, 1);
J(4, 2) = dGNdtheta2(4, 1);
clear bhat_gmm_theta2 dGNdtheta2;

% Obtaining the (1, 3), (2, 3), (3, 3), and (4, 3) components
bhat_gmm_theta3 = bhat_gmm;
bhat_gmm_theta3(3, 1) = bhat_gmm(3, 1)*1.001;
dGNdtheta3 = (GN(bhat_gmm_theta3) - GN(bhat_gmm))/(0.001*bhat_gmm(3, 1));
J(1, 3) = dGNdtheta3(1, 1);
J(2, 3) = dGNdtheta3(2, 1);
J(3, 3) = dGNdtheta3(3, 1);
J(4, 3) = dGNdtheta3(4, 1);
clear bhat_gmm_theta3 dGNdtheta3;
Gamma = abs(J);

% Obtaining GMM estimator var-cov matrix
varcov_matrix_gmm = inv(Gamma'*A_bhat_gmm*Gamma)*(Gamma'*A_bhat_gmm*inv(A_bhat_gmm)*A_bhat_gmm*Gamma)*inv(Gamma'*A_bhat_gmm*Gamma);
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
%% Questions asked in problem 3 that doesn't involve code
% Part 3a:
% PLEASE SEE THE SCANNED HANDWRITTEN PAGE FOR THE ANSWER..

% Part 3c:
% Attempting to estimate this model with only a mean independence
% assumption would be difficult due to what the integrating-out process
% involves (both complete and partial). In the context of me choosing to
% integrate out e1_i, in order to approximate the likelihood contribution
% of i (after integrating out e1_i), you need to know what the density of
% e1_i is up to parameters. If one does not, then one cannot take simulated
% draws of the error term being integrated out. Because of that, when
% estimating a non-linear and non-invertible model with ML, you need to
% assume the distribution of the error terms up to parameters(especially
% the one being integrated out).

% Part 3d:
% The computational results from part 3d code above shows that the
% difference between the first stage standard errors and second stage
% standard errors is negligible. This seems in indicate that the second
% stage process is more about finding the "correct" parameter estimates,
% and less about decreasing the standard errors of said estimates.

% Part 3f:
% The main change we see done to the model is that there is now an
% additional error term, eta_i, to the optimal price equation, p_i*. In
% terms of intuitive changes, because eta_i is observable to neither firm,
% they will essentially set the new error term to equal zero in their FOCs.
% This means when it comes to profit maximisation of either firm, the
% optimal quantity will not change. However, optimal p_i* is no longer
% pinned down by the optimal quantities of both firms due to the existence
% of error term eta. This could be interpreted as the optimal price being
% subject to potential demand shocks not accounted for by either firm. Note
% that this interpretation occurs when parameter phi is large. However, if
% phi is close to zero, the model is essentially identical to the original
% model.
% In terms of changes to estimation procedures, we would essentially be
% able to partially integrate out both e_1 and e_2. This means we would
% need to take simulated draws for said error terms when doing ML
% (2 sets of simulated draws instead of just 1 in the original model).
% Furthermore, we know have four parameters to estimate (theta1, theta2,
% sigma, and phi).
% For GMM, because we have four parameters to estimate now and we have four
% moments, our model can now be uniquely identified. However, we will need
% to do three sets of simulated draws since the model now has three error
% terms.

% Part 3g:
% With the introduction to new parameter rho and error term alpha in each
% firm's marginal cost functions, this means the FOCs and thus profit
% maximisations of both firms will change. As a result, the general optimal
% quantities chosen by both firms will be different from the original
% model. This is because the firms are now basically subject to potential
% supply-side shocks which affect their marginal costs. This is the case
% for cases when parameter rho is large. When rho is small, the model will
% behave closely to the original model.
% In terms of estimation changes, we now have four parameters to estimate
% in ML. Because we are able to partially integrate out both original error
% terms, we then need to take simulated draws for both error terms when
% doing ML.
% When doing GMM, we now have four parameters to estimate. Our modified
% model will still have four moments like the original model. Therefore,
% our model will be uniquely identified (like the changed model in part
% 3f). We again need to do three sets of simulated draws for error terms
% e_1, e_2, and alpha.

% Part 3h:
% All the following analyses are doing referring to the model in part 3f.
% Suppose we are in the scenario where we observe p_i, x1_i, x2_i, q1_i,
% and q2_i. Assume we see that q_i's are similar and x_i's are similar
% across markets. We can infer that sigma is small. If we further assume
% that the observed p_i's are quite different and fluctuating across
% markets, we can additionally conclude that phi will be large.
% similar, we can infer that sigma is small. Conversely, if we observe
% q_i's and x_i's to be different across markets by a large degree,
% this implies sigma is large. Additionally, if we see that p_i doesn't
% differ a lot, we can conclude that phi is small. Overall, we can say that
% observing this many variables in the data allows us to pin down phi and
% sigma from just looking at the data.
% Now assume we only observe total quantity q_i, p_i, x1_i, and x2_i.
% Assume we see that x1_i and x2_i are similar to each other across
% markets, and we also see that q_i's across markets are similar. Because
% we can't see the individual firm quantities, it is difficult to say from
% the data alone if sigma is small from our observations. However, if we
% see that p_i's are similar across markets, we can still infer that this
% is a world where phi is small (we can conclude the opposite if we see
% p_i's to differ by a lot across markets). Overall, this means that this
% change in what variables we can observe results in us to be able to pin
% down phi, but not sigma.
% Now assume we only observe p_i, x1_i, and x2_i (our original dataset). If
% we see that x_i's are similar across markets and p_i's are as well, this
% doesn't really tell us the value of sigma. This is because whilst x_i's
% are similar across markets, it's entirely possible for individual
% quantities to differ by a lot (but isn't observable to us). Furthermore,
% even if we see similar p_i's across markets, we can't pin down if this is
% due to being in a world with small sigma only, small phi only, or both.

% Part 3i:
% All the following analyses are doing referring to the model in part 3g.
% Assume we observe p_i, x1_i, x2_i, q1_i, and q2_i. Observe that both
% marginal cost equations share the same shock alpha and parameter rho.
% This means both individual quantities share the same alpha and rho as
% well. If we assume that the individual quantities and x_i's are similar
% across markets between firms, this means we can conclude sigma is small.
% However, we are unable to pin down the value/magnitude of rho.
% Conversely, if we see the individual quantities are different across
% markets (and the x_i's are different too), we can infer that sigma is
% large. But again, the value of rho cannot be pinned down. 
% Assume we observe p_i, x1_i, x2_i, and q_i. Now assume p_i's, x_i's, and
% q_i's are similar across markets. Due to only observing total quantity,
% we are unable to pin down the value of sigma. We also run into the issue
% of because alpha is shared in both individual marginal costs and
% quantities, we cannot pin down rho as well.
%==========================================================================