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
%==========================================================================

%==========================================================================
% Creating a function that calculates the likelihood contribution of
% consumer i, after integrating out \alpha_i, for each simulated draw of
% \alpha_i, and then sums them together along S. More specifically, this
% function calculates the inner probability of the product as seen in
% equation (3) of the problem set for customer i, then sums them together
% along S.
%=====
% NOTE
%=====
% WE ARE CHOOSING TO INTEGRATE OUT \alpha_i.

% Parameters:
% theta(1) = theta0;
% theta(2) = theta1;
% theta(3) = theta2;
% theta(4) = sigma_alpha
%=========
% END NOTE
%=========

function rhs_sum_component = l_io_alpha_i_i_sum(theta, i)
  global p_it y_it y_itm1;
  global S;
  global alpha;
  
  % When examining the RHS of equation (3) of the problem set, we see that
  % the probability inside of the product operator, for each consumer i and
  % each simulated draw of alpha_i, has a length of T. Applying the product
  % operator to this inner probability vector will collapse it to a scalar.
  % Initialising RHS sum of likelihood function.
  rhs_sum_component = 0;
  
  for s = 1:S
    % Creating U_it
    U_it = @(theta) exp(theta(1) + theta(2)*p_it(:, :, i) + theta(3)*y_itm1(:, :, i) + abs(theta(4))*alpha(:, s, i));
    % Creating probability inside of the product operator in the RHS of the
    % likelihood function.
    inner_p = @(theta) y_it(:, :, i).*(U_it(theta)./(1 + U_it(theta))) + (1 - y_it(:, :, i))./(1 + U_it(theta));
    % Applying product operator
    prod_inner_p = @(theta) prod(inner_p(theta));
    rhs_sum_component = rhs_sum_component + subsref(prod_inner_p(theta), struct('type', '()', 'subs', {{1, 1}}));
  end
  clear s;
end
%==========================================================================