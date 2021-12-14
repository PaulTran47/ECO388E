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

% ECO388E Problem Set 3; 3, 4
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
%==========================================================================

%==========================================================================
%% Part 3, 4: "Alternative-specific" value function iteration
%=====
% NOTE
%=====
% Like Rust (1987), we are using the conditional independence assumption.
% This means that \epsilon_{it} is noise that doesn't affect the future.
% Furthermore, recall that we assumed \epsilon_{it} are iid logit errors.
% As a result, we are able to analytically write out the expectation of the
% max in the following "alternative-specific" value functions.

% In terms of format, the overall value function is a 5x2 array, where the
% rows represent values of a_{t} and the columns represent values based on
% i_{t}. Because we have two "value" functions, \bar{v_{0}}(a_{t}) and
% \bar{v_{1}}(a_{t}), we have two columns.
%=========
% END NOTE
%=========
% Initialising discounting parameter and Euler's constant
beta = 0.9;
gamma = 0.5772;

% Creating vector that will house all possible values of a_{t}. Recall that
% once a machine reaches age 5, it will stay at age 5 forever until
% replaced.
a = (1:5)';

% Creating initial values for both value functions. For simplicity, we are
% making the value of not replacing or replacing, in all ages, to equal one
% initially.
v0 = ones(length(a), 1);
v1 = ones(length(a), 1);

%======================================================================
% Problem 4 code (setting parameter values to actually do the VFI loop)
%======================================================================
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

% Creating variable that stores max(abs(v - v_next)), the error metric. We
% are focusing on the max of the error because if the max value satisfies
% the threshold, every other value of the vector will.
error = 1;

% Creating variable that stores iteration number
iteration = 0;

%==============================
% Value function iteration loop
%==============================
% Iterating to find fixed point of value function. The threshold that stops
% the loop will be if the max error of the individual errors for v0 and v1
% is <= 0.001.
while error >= 0.001
  % Displaying iteration number
  disp(iteration);

  % Creating next-period age value function vectors that go into the RHS of
  % the value functions created by using initial or previous value
  % functions
  v0a1(1: length(a) - 1, 1) = v0(2:length(a), 1);
  v0a1(length(a), 1) = v0(length(a), 1);
  v1a1(1: length(a) - 1, 1) = v1(2:length(a), 1);
  v1a1(length(a), 1) = v1(length(a), 1);

  % Creating value functions created by using initial or previous value
  % functions
  %=====
  % NOTE
  %=====
  % Paramater value mapping:
  % theta(1) = \theta_{1}
  % theta(2) = R;
  %=========
  % END NOTE
  %=========
  v0_next = @(theta) theta(1).*a + beta.*(gamma + log(exp(v0a1) + exp(v1a1)));
  % Using repmat so the output of function v1_next is the same dimension as
  % vector of ages a.
  v1_next = @(theta) theta(2) + beta.*(gamma + log(exp(repmat(v0(1, 1), length(a), 1)) + exp(repmat(v1(1, 1), length(a), 1))));
  
  % Storing max(abs(v - v_next)). Choosing the max so the abs vector has
  % every value to be below error threshold.
  error = max(max(abs(v0 - v0_next(theta0))), max(abs(v1 - v1_next(theta0))));
  disp(error);
  
  % Storing previous value function
  v0 = v0_next(theta0);
  v1 = v1_next(theta0);
  
  iteration = iteration + 1;
end
clear iteration error v0a1 v1a1;
%==========================================================================

%==========================================================================
%% Part 4: Using VFI loop results when given parameter values \theta_{1} = -1, R = -3 to answer questions
%=====
% NOTE
%=====
% We focus on age a_{t} = 2. To see whether or not a firm replaces its
% machine, we need to calculate \epsilon_{0t} - \epsilon_{1t}. This is
% equivalent to calculating the difference between our "value functions":
% \bar{v_{0}}(a_{t}) - \bar{v_{1}}(a_{t}). If a firm is indifferent with
% replacing or keeping a machine, this means

% \bar{v_{0}}(2) = \bar{v_{1}}(2)
% \implies \theta_{1}*2 + \epsilon_{0t} + \beta*E[v(3, \epsilon_{3}; \theta_{t + 1})] = R + \epsilon_{0t} + \beta*E[v(1, \epsilon_{1}; \theta_{t + 1})]
% \implies \epsilon_{0t} - \epsilon_{1t} = \theta_{1}*2 - R + \beta*E[v(3, \epsilon_{3}; \theta_{t + 1}) - v(1, \epsilon_{1}; \theta_{t + 1})].

%=========
% END NOTE
%=========
% Calculating \epsilon_{0t} - \epsilon_{1t}
a = 2;
diff_e0e1 = theta0(1, 1)*a - theta0(1, 2) + beta*((gamma + log(exp(v0(a + 1, 1)) + exp(v1(a + 1, 1)))) - (gamma + log(exp(v0(1, 1)) + exp(v1(1, 1)))));

%=======
% ANSWER
%=======
% Given our calculated difference in the error terms, we see that a firm is
% indifferent between replacing and keeping a machine whose age is 2 when
% said difference is equal to -0.114487971844367. If the difference is less
% than -0.114487971844367, then the firm will not replace the machine. If
% the difference is greater than -0.114487971844367, then the firm will
% replace the machine.
%===========
% END ANSWER
%===========

% Calculating the probability that the firm replaces its machine (using the
% extreme value distribution given our assumptions about \epsilon_{it}
% being iid logit errors)
p1_diff_e0e1 = exp(-diff_e0e1)/(1 + exp(-diff_e0e1));

%=======
% ANSWER
%=======
% Using the extreme value distribution to calculate the probability of
% seeing this difference (which results in seeing the firm replace its
% machine at age 2), the probability is 0.528590770331370.
%===========
% END ANSWER
%===========

% Calculating PDV of future profits for a firm at state {a_{t} = 4,
% \epsilon_{0t} = 1, \epsilon_{1t} = -1.5}
% Initialising specific variables
a = 4;
e0 = 1;
e1 = -1.5;
pdv = max(theta0(1, 1)*a + e0 + beta.*v0(a, 1), theta0(1, 2) + e1 + beta.*v1(a, 1));

%=======
% ANSWER
%=======
% The PDV of future profits for a firm at state {a_{t} = 4,
% \epsilon_{0t} = 1, \epsilon_{1t} = -1.5} is -14.7547873524102. As a
% result, we see that it is still cheaper for the firm to replace the
% machine in this period than to wait until age 5 to replace.
%===========
% END ANSWER
%===========
%==========================================================================