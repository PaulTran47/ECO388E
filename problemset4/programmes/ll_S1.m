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
% Creating a function that calculates the log of the kernel density
% estimator value when given a bandwidth h and individual S1 value,
% performs this function for every single S1 value, then add them all
% together. This sum will be the objective function for MLE.
function S1_sum = ll_S1(h)
  global S1 n1;
  global ll_S1_i;
  S1_sum = 0;
  for i = 1:n1
    S1_sum = S1_sum + ll_S1_i(S1(i), h);
  end
  clear i;
end
%==========================================================================