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
% Creating a function that creates the sum needed for simulating the
% inner expectation E[p|x1_i, x2_i, e1_i, e2_i; theta1, theta2, sigma] of
% the generic moments condition.
function rhs_sum_component = inner_expectations_p_sum(theta)
  global x1_i x2_i;
  global N S;
  global e_i_matrix;
  rhs_sum_component = zeros(N, 1);
  for i = 1:S
    % Choosing a simulated draw of e1_i
    e1_i = e_i_matrix(:, i);
    % Choosing a simulated draw of e2_i
    e2_i = e_i_matrix(:, 20 + i);
  
    % theta(3) = sigma
    % Creating the RHS of our model, which the function will use to
    % simulate the inner expectation mentioned above
    model_rhs = @(theta) (100 + exp(theta(1) + theta(2)*x1_i + abs(theta(3))*e1_i) + exp(theta(1) + theta(2)*x2_i + abs(theta(3))*e2_i))/3;
  
    % Calculating integrand for the chosen simulated draw of e1_i, then
    % adding it to rhs_sum_component
    rhs_sum_component = rhs_sum_component + subsref(model_rhs(theta), struct('type', '()', 'subs', {{1:N, 1}}));
  end
  clear i;
end
%==========================================================================