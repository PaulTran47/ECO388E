%==========================================================================
%% 2 percentage signs represent comments for code;
% 1 percentage sign represents commented out code;

% Creating a function that calculates the likelihood contribution of i,
% after integrating out e1_i, for each simulated draw of e1_i, and then
% sums them together along S.
%=======================================
% WE ARE CHOOSING TO INTEGRATE OUT e1_i:
%=======================================
function rhs_sum_component = l_io_e1_i_i_sum(theta)
  global x1_i x2_i regressand;
  global N S;
  global e1_i_matrix;
  rhs_sum_component = zeros(N, 1);
  for i = 1:S
    % Choosing a simulated draw of e1_i
    e1_i = e1_i_matrix(:, i);
  
    % theta(3) = sigma
    % Creating inverse function e2_i
    % We are also making sure there are no negative values by replacing
    % then with 0.0001.
    e2_i = @(theta) (log(max(3*regressand - 100 - exp(theta(1) + theta(2)*x1_i + abs(theta(3))*e1_i), 0.0001)) - theta(1) - theta(2)*x2_i)/abs(theta(3));
    % Creating the likelihood contribution of i, after integrating out e1_i
    l_io_e1_i_i = @(theta) normpdf(e2_i(theta)).*abs(3./(abs(theta(3))*(3*regressand - 100 - exp(theta(1) + theta(2)*x1_i + abs(theta(3))*e1_i))));
  
    % Calculating l_io_e1_i_i for the chosen simulated draw of e1_i, then
    % adding it to rhs_sum_component
    rhs_sum_component = rhs_sum_component + subsref(l_io_e1_i_i(theta), struct('type', '()', 'subs', {{1:N, 1}}));
  end
  clear i;
end
%==========================================================================