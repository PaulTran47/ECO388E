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
% Creating a function that will calculate the MSE between the individual
% true wages in the testing set, for every combination of education and
% gender found in the testing set.
function mse_sum = mse(lambda)
  global test test_unique
  % Initialising scalar holding MSE sum
  mse_sum = 0;
  for i = 1:length(test_unique)
    for j = [0, 1]
      mse_sum_component = sum((test(test(:, 1) == i, 2) - subsref(knn(lambda), struct('type', '()', 'subs', {{i, j + 1}}))).^2);
      mse_sum = mse_sum + mse_sum_component;
    end
  end
  mse_sum = mse_sum/length(mse_sum);
end
%==========================================================================