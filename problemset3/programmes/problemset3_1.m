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

% ECO388E Problem Set 3, 1
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
%% Part 1: Write down the dynamic programming problem for a firm maximising the PDF of future profits (assume an infinite horizon)
%=======
% ANSWER
%=======
% Observe that the state variables needed to define the current period's
% profits and expectations about discounted future profits are
% {a_{t}, \epsilon_{t}}. This means we can write the firm's value function
% as follows:

% V(a_{t}, \epsilon_{t}; \theta) = max_{i_{t}}(E[\sum_{j = t}^{\infty} \beta^{j - t}\Pi(a_{j}, i_{j}, \epsilon_{j},; \theta) | a_{t}, i_{t}; \theta]).

% Converting this into a Bellman equation yields the following:

% V(a_{t}, \epsilon_{0t}, \epsilon_{1t}; \theta) = max_{i_{t}}(\Pi(a_{t}, \epsilon_{0t}, \epsilon_{1t}; \theta) + \beta E[V(a_{t + 1}, \epsilon_{0t + 1}, \epsilon_{1t + 1}; \theta) | a_{t}, i_{t}; \theta]).
%===========
% END ANSWER
%===========
%==========================================================================