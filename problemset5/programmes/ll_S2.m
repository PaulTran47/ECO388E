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
% estimator value when given a bandwidth h and individual S2 value,
% performs this function for every single S2 value, then add them all
% together. This sum will be the objective function for MLE.
function S2_sum = ll_S2(h)
  global S2 n2;
  global ll_S2_i;
  S2_sum = 0;
  for i = 1:n2
    S2_sum = S2_sum + ll_S2_i(S2(i), h);
  end
  clear i;
end
%==========================================================================