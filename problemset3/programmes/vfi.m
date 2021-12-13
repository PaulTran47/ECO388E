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
% Creating a function that performs value function iteration (VFI) for Rust
% (1987)'s "alternative-specific" value functions.
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
function v0v1_matrix = vfi(theta)
  % Calling global variables to be used
  global beta gamma a_max;

  % Creating vector that will house all possible values of a_{t}. Recall
  % that once a machine reaches age 5, it will stay at age 5 forever until
  % replaced. This vector essentially represents the allowed values age can
  % be.
  a_vfi = (1:5)';

  % Creating initial values for both value functions. For simplicity, we
  % are making the value of not replacing or replacing, in all ages, to
  % equal one initially.
  v0 = ones(a_max, 1);
  v1 = ones(a_max, 1);
  
  % Creating variable that stores max(abs(v - v_next)), the error metric.
  % We are focusing on the max of the error because if the max value
  % satisfies the threshold, every other value of the vector will.
  error = 1;

  % Creating variable that stores iteration number
  iteration = 0;
  
  %==============================
  % Value function iteration loop
  %==============================
  % Iterating to find fixed point of value function. The threshold that
  % stops the loop will be if the max error of the individual errors for v0
  % and v1 is <= 0.001.
  while error >= 0.001
    % Displaying iteration number
    disp(iteration);

    % Creating next-period age value function vectors that go into the RHS
    % of the value functions created by using initial or previous value
    % functions
    v0a1(1: a_max - 1, 1) = v0(2:a_max, 1);
    v0a1(a_max, 1) = v0(a_max, 1);
    v1a1(1: a_max - 1, 1) = v1(2:a_max, 1);
    v1a1(a_max, 1) = v1(a_max, 1);

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
    v0_next = @(theta) theta(1).*a_vfi + beta.*(gamma + log(exp(v0a1) + exp(v1a1)));
    % Using repmat so the output of function v1_next is the same dimension
    % as vector of ages a.
    v1_next = @(theta) theta(2) + beta.*(gamma + log(exp(repmat(v0(1, 1), a_max, 1)) + exp(repmat(v1(1, 1), a_max, 1))));
  
    % Storing max(abs(v - v_next)). Choosing the max so the abs vector has
    % every value to be below error threshold.
    error = max(max(abs(v0 - v0_next(theta))), max(abs(v1 - v1_next(theta))));
    disp(error);
  
    % Storing previous value function
    v0 = v0_next(theta);
    v1 = v1_next(theta);
  
    % Creating matrix that stores v0 and v1 as columns in that order
    v0v1_matrix = [v0 v1];
  
    iteration = iteration + 1;
  end
  clear iteration error v0a1 v1a1;
end
%==========================================================================