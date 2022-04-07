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
% Creating a function performs KNN search, given K (function of \lambda)
% for every unique education value in the test subset. Specifically, for
% each education value, KNN search is performed using the training subset.
% Then, the mean will be performed on the corresponding wage values in the
% training subset. This will finally be placed into a vector, which is the
% final output.
function knn_predictions = knn(lambda)
  global train;
  %global test;
  global test_unique;
  global K;

  % Declaring vector that will hold associated estimate of
  % E[wage | education, gender].
  m_knn_test = zeros(length(test_unique), 2);
  
  for i = 1:length(test_unique)
    for j = [0, 1]
      %disp(test(i, 1));
      Idx = knnsearch(train(:, 1:2), [test_unique(i) j], 'K', round(K(lambda)));
    
      % Finding average of corresponding wage values to these neighbours
      % (from the training subset)
      wage_train = train(:, 2);
      m_knn_test(i, j + 1) = mean(wage_train(Idx, :));
    end
  end
  
  knn_predictions = m_knn_test;
end
%==========================================================================