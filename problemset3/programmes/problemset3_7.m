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

% ECO388E Problem Set 3, 7
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
%% Part 7a: Two types of firms differing in their value of \theta_{1}.
%=======
% ANSWER
%=======
% We now assume that two types of firms exist, where they differ in their
% value of \theta_{1}. Furthermore, we assume that proportion \alpha of
% firms have \theta_{1} = \theta_{1A} and proportion (1 - \alpha) have
% \theta_{1} = \theta_{1B}. In order to accommodate these changes, we need
% to modify our dynamic programming problem -- specifically, the Bellman
% equation -- as follows:

% V_{j}(a_{jt}, \epsilon_{0jt}, \epsilon_{1jt}; \theta_{j}) = max_{i_{jt}}(\Pi(a_{jt}, \epsilon_{0jt}, \epsilon_{1jt}; \theta_{j}) + \beta E[V(a_{jt + 1}, \epsilon_{0jt + 1}, \epsilon_{1jt + 1}; \theta_{j}) | a_{jt}, i_{jt}; \theta_{j}]).

% Observe that the value function, replacement choice shocks, replacement
% decision, machine age, and parameters now have a subscript j. This is to
% represent that two types of firms exist now. Essentially, this means that
% our dynamic programming problem now involves calculating two different
% value functions via value function iteration.

% Changes that need to be applied to the likelihood function are as
% follows:

% p(i_{j} = 1 | a_{jt}, \epsilon_{jt}) = \alpha*p(i_{j} = 1 | a_{jt}, \epsilon_{jt}; \theta_{1A}) + (1 - \alpha)*p(i_{j} = 1 | a_{jt}, \epsilon_{jt}; \theta_{1B}).

% If we assume both the error terms are iid logit errors and conditional
% independence, we obtain the expression for the likelihood function:

% p(i_{j} = 1 | a_{jt}, \epsilon_{jt}) = \alpha*\frac{e^{V_{1A}(a_{jt}, \epsilon_{jt})}}{e^{V_{1A}(a_{jt}, \epsilon_{jt})} + e^{V_{0A}(a_{jt}, \epsilon_{jt})}} + (1 - \alpha)*\frac{e^{V_{1B}(a_{jt}, \epsilon_{jt})}}{e^{V_{1B}(a_{jt}, \epsilon_{jt})} + e^{V_{0B}(a_{jt}, \epsilon_{jt})}}.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 7b: Model described in part 7a, but we have panel data. Write down the likelihood function

%=======
% ANSWER
%=======
% Given our model written in part 7a, if we had panel data that ends at
% time period T, we must now consider the sequence of decisions by a firm,
% {i(a_{jt}, \epsilon_{jt}) = 1}_{t < T}. From this, we can
% write the likelihood function as follows:

% p({i(a_{jt}, \epsilon_{jt}) = 1}_{t < T}) = \alpha*p({i(a_{jt}, \epsilon_{jt}) = 1}_{t < T}; \theta_{1A}) + (1 - \alpha)*p({i(a_{jt}, \epsilon_{jt}) = 1}_{t < T}; \theta_{1B}).

% As before, if we assume both the error terms are iid logit errors and
% conditional independence, we obtain the following expression for the
% likelihood function:

% p({i(a_{jt}, \epsilon_{jt}) = 1}_{t < T}) = \alpha*\prod_{t < T} \frac{e^{V_{1A}(a_{jt}, \epsilon_{jt})}}{e^{V_{1A}(a_{jt}, \epsilon_{jt})} + e^{V_{0A}(a_{jt}, \epsilon_{jt})}} + (1 - \alpha)*\prod_{t < T} \frac{e^{V_{1B}(a_{jt}, \epsilon_{jt})}}{e^{V_{1B}(a_{jt}, \epsilon_{jt})} + e^{V_{0B}(a_{jt}, \epsilon_{jt})}}.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 7c: Given assumes made in parts 7a and 7b, assume machines differ in \theta_{1}.
%=======
% ANSWER
%=======
% Given all the assumptions made in parts 7a and 7b, if we additionally
% assume that machines now can differ in parameter \theta_{1} in the same
% split as firms do (i.e., a firm's new machine could have \theta_{1A} with
% probability \alpha or \theta_{1B} with probability (1 - \alpha)), this is
% essentially saying that the model now considers machine type to be a
% state variable that has to be included. Additionally,, machines differing
% will now result in our error terms \epsilon_{jt} to be serially
% correlated. This will result in the expectations operator seen in all our
% "alternative-specific" value functions to be conditioned on \epsilon_{jt}
% now. Due to there being many "alternative-specific" value functions now
% with the different types assumed, we will write the "general" change that
% all these value functions experience. More specifically, we have:

% \bar{V}_{0j} = u(a_{jt}, 0; \theta_{j}) + \beta*E[V(a_{jt + 1}, \epsilon_{jt + 1}; \theta_{j}) + a_{jt}, \epsilon_{jt}, i_{jt} = 0; \theta_{j}];
% \bar{V}_{1j} = u(a_{jt}, 1; \theta_{j}) + \beta*E[V(a_{jt + 1}, \epsilon_{jt + 1}; \theta_{j}) + a_{jt}, \epsilon_{jt}, i_{jt} = 1; \theta_{j}].

% With regards to changes in our likelihood calculation, recall that
% machines differing in type translates to the choice-specific error terms
% to now be serially correlated to each other. As a result, the formula
% used to calculate the likelihood (as done in problems 3-6 and parts 7a
% and 7b) is no longer applicable. In terms of class notes, this is because
% assumption C is violated.
%===========
% END ANSWER
%===========
%==========================================================================

%==========================================================================
%% Part 7d: Initial conditions problem
%=======
% ANSWER
%=======
% In simple terms, the initial conditions problem that was ignored in parts
% 7a - 7c essentially comes from the fact that we do not see the initial
% period of data. In other words, whilst we have a "complete" model of how
% (i_{j1}, \ldots, i_{jT}) and (a_{j2}, \ldots, a_{jT}) are determined, we
% do not have a "complete" model of how age a_{j1} is determined. This is
% because we do not observe a_{j0}, \epsilon_{j0}, and i_{0t}. As a result,
% any serial correlation in the error terms that our model observes will
% either be due to unobserved heterogeneity or by the initial draw of
% \epsilon_{jt}. We should note that this is a problem just for the
% econometrician.

% One solution to the initial conditions problem is to condition on a_{j1},
% which yields the following likelihood:

% p(a_{j2}, \ldots, a_{jT}, i_{j1}, \ldots, i_{jT} | a_{j1}; \theta_{j}).

% Doing this will require us to know the density

% p(\epsilon_{1} | a_{1}; \theta).

% The solution would then require simulating the above conditional
% distribution given the solution of the dynamic programming problem.

% An additional potential solution could be to simulate \epsilon_{jt} all
% the way back to the initial period of data.
%===========
% END ANSWER
%===========

%==========================================================================

%==========================================================================
%% Part 7e: Evolution of age is stochastic
%=======
% ANSWER
%=======
% Assume the evolution of a_{t} is now the following if you don't replace:

% a_{t + 1} = 
%  min(5, a_{t} + 1), with probability \lambda;
%  a_{t}, with probability (1 - \lambda).

% It should be noted that if a firm does replace, the machine's age will
% still always be one. This all means that a machine's age isn't guaranteed
% to increase by one if it is kept by the firm into the next period.
% We also assume that our data is a random sample of firms that have
% existed for a long time. Despite the stochastic nature of a machine's age
% if it isn't replaced, there is still some information in the data that
% can be used to figure out \lambda. Essentially, if we see firms are
% waiting longer than expected to replace their machines via the machines'
% ages seen in the data, we can use this information to determine \lambda.
% In particular, by using the distribution of the unobservables,
% \epsilon_{it}, and other existing assumptions and information made and
% given about our original model, we can calculate the expected age of when
% a firm should replace its machine. Therefore, if we see across all
% observations/firms that engines have an age older than the expected value
% (which is when an engine is expected to be replaced), this will allow us
% to learn \lambda.
%===========
% END ANSWER
%===========
%==========================================================================