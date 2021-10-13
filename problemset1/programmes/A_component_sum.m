%==========================================================================
%% 2 percentage signs represent sections of code;
% 1 percentage sign represents comments for code or commented out code;

% Creating a function that find the sum needed for the estimate of the
% optimal weight matrix A. How it works is that it takes the individual
% component from the component matrix made above, multiplies it with its
% transpose, then finds the sum for all i.
function component_i_component_i_t_sum = A_component_sum(theta)
  global N;
  global component;
  % Storing the needed dimension of component
  A_dim = min(size(component(theta)));
  component_i_component_i_t_sum = zeros(A_dim, A_dim);
  for i = 1:N
    component_i = subsref(component(theta), struct('type', '()', 'subs', {{1:A_dim, i}}));
    component_i_component_i_t = component_i*component_i';
    component_i_component_i_t_sum = component_i_component_i_t_sum + component_i_component_i_t;
  end
  clear i;
end
%==========================================================================