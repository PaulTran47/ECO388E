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

% ECO388E Problem Set 3, 2
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
% \Pi(a_{t}, i_{t}, \epsilon_{0t}, \epsilon_{1t}) = 
%   \theta_{1}a_{t} + \epsilon_{0t}, if i_{t} = 0;
%   R + \epsilon_{1t}, if i_{t} = 1,
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
%% Part 2: What are the differences between this \Pi(a_{t}, i_{t}, \epsilon_{0t}, \epsilon_{1t}) and the profit function in the class notes on Rust? What happened to c(0; \theta)?
%=======
% ANSWER
%=======
% The first major difference between the problem set and Rust (1987)'s
% original problem (the HZ problem) is that in the latter, the bus engine
% mileage, x_{t}, is a continuous variable. In the former, the equivalent
% machine age, a_{t}, is a discrete variable that can only take values in
% the set {1, 2, 3, 4, 5}. The second difference we see between the profit
% functions of the problem set and the class notes is that the expected
% cost function of operating a machine with (a_{t}, c(a_{t}; \theta), is
% equal to \theta_{1}a_{t} in the problem set. As a result, we have in the
% problem set that c(0; \theta) = 0. This makes sense because having an age
% of a_{t} = 0 is not technically possible.
% Other technical differences are that the signs of R and c(a_{t}; \theta)
% are reversed between the two profit functions. However, this should only
% be notational differences, not actual-value differences. This is because
% both profit functions still treat the replacement cost, R, and expected
% cost function, c(.), to be negative.
%===========
% END ANSWER
%===========
%==========================================================================