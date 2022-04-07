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
% ECO388E Problem Set 5 (2, Spring 2022), 1
% Paul Le Tran, plt377
% 2 April, 2022
%==========================================================================

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'path\to\programmes';
data_dir = 'path\to\data';

cd(home_dir);
%==========================================================================

%==========================================================================
%% Importing data
cd(data_dir);
opts = delimitedTextImportOptions("NumVariables", 25);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["VarName1", "wage", "educ", "exper", "tenure", "nonwhite", "female", "married", "numdep", "smsa", "northcen", "south", "west", "construc", "ndurman", "trcommpu", "trade", "services", "profserv", "profocc", "clerocc", "servocc", "lwage", "expersq", "tenursq"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
wage1 = readtable(append(data_dir, '\wage1.csv'), opts);
% Clear temporary variables
clear opts

% Returning to home directory
cd(home_dir);
%==========================================================================

%==========================================================================
%% Part 1i: Use the kernel method to estimate the conditional expectation of wage on education and gender
%=====
% NOTE
%=====
% Because we are unsure of the problem's wording, we will for now estimate
% E[wage | education, gender].
%=========
% END NOTE
%=========
%====================================================================
% We are choosing bandwidth h based on the rule-of-thumb method. More
% information can be found at:
% https://en.wikipedia.org/wiki/Kernel_density_estimation#A_rule-of-thumb_bandwidth_estimator
% More specifically, we are using Silverman's rule-of-thumb.
% Creating vector of wages from data table
%====================================================================
wage = table2array(wage1(:, 2));
% Creating vector of education from data table
edu = table2array(wage1(:, 3));
% Creating vector of gender (female = 1, male = 0) from data table
gender = table2array(wage1(:, 7));

% Getting sample size
n = length(wage);
% Calculating standard deviation of sample wages
sigma_hat = std(edu);
% Calculating bandwidth h according to Silverman's rule-of-thumb and with
% respect to education
h = 0.9*min(sigma_hat, (iqr(edu)/1.34))*(n^(-1/5));

% Creating Gaussian kernel function (input is education)
K_edu = @(x_edu) normpdf((x_edu - edu)./h);

% Creating indicator function for female (no kernel function because gender
% is binary/discrete)
indicator_female = gender == 1;

% Creating indicator function for male (no kernel function because gender
% is binary/discrete)
indicator_male = gender == 0;

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are education and gender (female))
m_numerator_female = @(x_edu) sum(wage.*K_edu(x_edu).*indicator_female);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are education and gender (female))
m_denominator_female = @(x_edu) sum(K_edu(x_edu).*indicator_female);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% education and gender(female))
m_female = @(x_edu) m_numerator_female(x_edu)./m_denominator_female(x_edu);

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are education and gender (male))
m_numerator_male = @(x_edu) sum(wage.*K_edu(x_edu).*indicator_male);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are education and gender (male))
m_denominator_male = @(x_edu) sum(K_edu(x_edu).*indicator_male);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% education and gender(female))
m_male = @(x_edu) m_numerator_male(x_edu)./m_denominator_male(x_edu);

% Creating kernel conditional expectation estimation plots
fplot(m_female, [min(edu) max(edu)], 'LineWidth', 2)
hold on
fplot(m_male, [min(edu) max(edu)], 'LineWidth', 2)
grid on
xlabel('Education');
ylabel('m(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = 'gender = female = 1';
Legend{2} = 'gender = male = 0';
legend(Legend, 'Location', 'best');

hold on
title({'Kernel conditional expectation estimator', '$\hat{m}_{h}(x) \approx E[wage | education, gender]$,', 'Using Gaussian kernel,' ,'h calculated using the Silverman rule-of-thumb.'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\1ia_plot.png');
close(gcf);

%=======================================================
% We are choosing bandwidth h based on cross-validation.
%=======================================================
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Creating vector of education from data table
edu = table2array(wage1(:, 3));
% Creating vector of education from data table. For cross-validation
% purposes, we are randomising the education vector to obtain better
% distributed testing and training sets.
edu_rand = edu(randperm(length(edu)));
% Creating vector of gender (female = 1, male = 0) from data table
gender = table2array(wage1(:, 7));

% Getting sample size
n = length(edu_rand);
% Creating subsample S1
global S1;
S1 = edu_rand(1:(n/2), 1);
% Getting S1 size
global n1;
n1 = length(S1);
% Creating subsample S2
global S2;
S2 = edu_rand(((n/2) + 1):n, 1);
% Getting S2 size
global n2;
n2 = length(S2);

%=====
% NOTE
%=====
% Obtaining MLE of S2 sample
%=========
% END NOTE
%=========
% Creating kernel function that takes a single S2 value as input
global K_S2_i;
K_S2_i = @(x, h) normpdf((S1 - x)./h);

% Creating density function that takes a single S2 value as input
global f_S2_i;
f_S2_i = @(x, h) sum(K_S2_i(x, h))./(n1.*h);

% Creating log-likelihood contribution of a single S2 value using above
% density function
global ll_S2_i;
ll_S2_i = @(x, h) log(f_S2_i(x, h));

% Creating objective function for MLE (takes S2 values as input for
% calculation). More information about the ll_S2 function can be found in
% ll_S2.m
objfcn_S2_ml = @(h) ll_S2(h);

% Initial bandwidth value for numerical optimisation (setting equal to
% Silverman rule-of-thumb)
h0 = 0.383681233738354;

% Numerical optimisation to find bandwidth h (using S2 values as input)
h_S2_ml = fminunc(@(h) -objfcn_S2_ml(h), h0);

%=====
% NOTE
%=====
% Obtaining MLE of S1 sample
%=========
% END NOTE
%=========
% Creating kernel function that takes a single S1 value as input
global K_S1_i;
K_S1_i = @(x, h) normpdf((S2 - x)./h);

% Creating density function that takes a single S1 value as input
global f_S1_i;
f_S1_i = @(x, h) sum(K_S1_i(x, h))./(n2.*h);

% Creating log-likelihood contribution of a single S1 value using above
% density function
global ll_S1_i;
ll_S1_i = @(x, h) log(f_S1_i(x, h));

% Creating objective function for MLE (takes S1 values as input for
% calculation). More information about the ll_S1 function can be found in
% ll_S1.m
objfcn_S1_ml = @(h) ll_S1(h);

% Initial bandwidth value for numerical optimisation (setting equal to the
% Silverman rule-of-thumb value).
h0 = 0.383681233738354;

% Numerical optimisation to find bandwidth h (using S1 values as input)
h_S1_ml = fminunc(@(h) -objfcn_S1_ml(h), h0);

%=====
% NOTE
%=====
% Currently using the average of h_S1_ml and h_S2_ml.
% Because of the random permutation applied on the wage vector at the
% beginning, MLE estimations using S1 and S2 subsamples always change. The
% reference estimates used for the following density plot are as follows:
% h_S1_ml = 0.753840599172369;
% h_S2_ml = 0.184990600751972.
%=========
% END NOTE
%=========
h = (h_S1_ml + h_S2_ml)/2;

% Creating Gaussian kernel function (input is education)
K_edu = @(x_edu) normpdf((x_edu - edu)./h);

% Creating indicator function for female (no kernel function because gender
% is binary/discrete)
indicator_female = gender == 1;

% Creating indicator function for male (no kernel function because gender
% is binary/discrete)
indicator_male = gender == 0;

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are education and gender (female))
m_numerator_female = @(x_edu) sum(wage.*K_edu(x_edu).*indicator_female);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are education and gender (female))
m_denominator_female = @(x_edu) sum(K_edu(x_edu).*indicator_female);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% education and gender(female))
m_female = @(x_edu) m_numerator_female(x_edu)./m_denominator_female(x_edu);

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are education and gender (male))
m_numerator_male = @(x_edu) sum(wage.*K_edu(x_edu).*indicator_male);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are education and gender (male))
m_denominator_male = @(x_edu) sum(K_edu(x_edu).*indicator_male);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% education and gender(female))
m_male = @(x_edu) m_numerator_male(x_edu)./m_denominator_male(x_edu);

% Creating kernel conditional expectation estimation plots
fplot(m_female, [min(edu) max(edu)], 'LineWidth', 2)
hold on
fplot(m_male, [min(edu) max(edu)], 'LineWidth', 2)
grid on
xlabel('Education');
ylabel('m(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = 'gender = female = 1';
Legend{2} = 'gender = male = 0';
legend(Legend, 'Location', 'best');

hold on
title({'Kernel conditional expectation estimator', '$\hat{m}_{h}(x) \approx E[wage | education, gender]$,', 'Using Gaussian kernel,' ,'h is the average of MLE estimates found', 'from cross-validation via splitting the sample', 'into equal-sized sub-samples S1 and S2'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\1ib_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 1ii: Use the kernel method to estimate the conditional expectation of wage on experience and gender
%=====
% NOTE
%=====
% Because we are unsure of the problem's wording, we will for now estimate
% E[wage | experience, gender].
%=========
% END NOTE
%=========
%====================================================================
% We are choosing bandwidth h based on the rule-of-thumb method. More
% information can be found at:
% https://en.wikipedia.org/wiki/Kernel_density_estimation#A_rule-of-thumb_bandwidth_estimator
% More specifically, we are using Silverman's rule-of-thumb.
%====================================================================
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Creating vector of experience from data table
exper = table2array(wage1(:, 4));
% Creating vector of gender (female = 1, male = 0) from data table
gender = table2array(wage1(:, 7));

% Getting sample size
n = length(wage);
% Calculating standard deviation of sample wages
sigma_hat = std(exper);
% Calculating bandwidth h according to Silverman's rule-of-thumb and with
% respect to experience
h = 0.9*min(sigma_hat, (iqr(exper)/1.34))*(n^(-1/5));

% Creating Gaussian kernel function (input is experience)
K_exper = @(x_exper) normpdf((x_exper - exper)./h);

% Creating indicator function for female (no kernel function because gender
% is binary/discrete)
indicator_female = gender == 1;

% Creating indicator function for male (no kernel function because gender
% is binary/discrete)
indicator_male = gender == 0;

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are experience and gender (female))
m_numerator_female = @(x_exper) sum(wage.*K_exper(x_exper).*indicator_female);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are experience and gender (female))
m_denominator_female = @(x_exper) sum(K_exper(x_exper).*indicator_female);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% experience and gender(female))
m_female = @(x_exper) m_numerator_female(x_exper)./m_denominator_female(x_exper);

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are experience and gender (male))
m_numerator_male = @(x_exper) sum(wage.*K_exper(x_exper).*indicator_male);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are experience and gender (male))
m_denominator_male = @(x_exper) sum(K_exper(x_exper).*indicator_male);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% experience and gender(female))
m_male = @(x_exper) m_numerator_male(x_exper)./m_denominator_male(x_exper);

% Creating kernel conditional expectation estimation plots
fplot(m_female, [min(exper) max(exper)], 'LineWidth', 2)
hold on
fplot(m_male, [min(exper) max(exper)], 'LineWidth', 2)
grid on
xlabel('Experience');
ylabel('m(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = 'gender = female = 1';
Legend{2} = 'gender = male = 0';
legend(Legend, 'Location', 'Northeast');

hold on
title({'Kernel conditional expectation estimator', '$\hat{m}_{h}(x) \approx E[wage | experience, gender]$,', 'Using Gaussian kernel,' ,'h calculated using the Silverman rule-of-thumb.'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\1iia_plot.png');
close(gcf);

%=======================================================
% We are choosing bandwidth h based on cross-validation.
%=======================================================
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Creating vector of experience from data table
exper = table2array(wage1(:, 4));
% Creating vector of experience from data table. For cross-validation
% purposes, we are randomising the experience vector to obtain better
% distributed testing and training sets.
exper_rand = exper(randperm(length(exper)));
% Creating vector of gender (female = 1, male = 0) from data table
gender = table2array(wage1(:, 7));

% Getting sample size
n = length(exper_rand);
% Creating subsample S1
global S1;
S1 = exper_rand(1:(n/2), 1);
% Getting S1 size
global n1;
n1 = length(S1);
% Creating subsample S2
global S2;
S2 = exper_rand(((n/2) + 1):n, 1);
% Getting S2 size
global n2;
n2 = length(S2);

%=====
% NOTE
%=====
% Obtaining MLE of S2 sample
%=========
% END NOTE
%=========
% Creating kernel function that takes a single S2 value as input
global K_S2_i;
K_S2_i = @(x, h) normpdf((S1 - x)./h);

% Creating density function that takes a single S2 value as input
global f_S2_i;
f_S2_i = @(x, h) sum(K_S2_i(x, h))./(n1.*h);

% Creating log-likelihood contribution of a single S2 value using above
% density function
global ll_S2_i;
ll_S2_i = @(x, h) log(f_S2_i(x, h));

% Creating objective function for MLE (takes S2 values as input for
% calculation). More information about the ll_S2 function can be found in
% ll_S2.m
objfcn_S2_ml = @(h) ll_S2(h);

% Initial bandwidth value for numerical optimisation (setting equal to
% Silverman rule-of-thumb)
h0 = 3.488946569653070;

% Numerical optimisation to find bandwidth h (using S2 values as input)
h_S2_ml = fminunc(@(h) -objfcn_S2_ml(h), h0);

%=====
% NOTE
%=====
% Obtaining MLE of S1 sample
%=========
% END NOTE
%=========
% Creating kernel function that takes a single S1 value as input
global K_S1_i;
K_S1_i = @(x, h) normpdf((S2 - x)./h);

% Creating density function that takes a single S1 value as input
global f_S1_i;
f_S1_i = @(x, h) sum(K_S1_i(x, h))./(n2.*h);

% Creating log-likelihood contribution of a single S1 value using above
% density function
global ll_S1_i;
ll_S1_i = @(x, h) log(f_S1_i(x, h));

% Creating objective function for MLE (takes S1 values as input for
% calculation). More information about the ll_S1 function can be found in
% ll_S1.m
objfcn_S1_ml = @(h) ll_S1(h);

% Initial bandwidth value for numerical optimisation (setting equal to the
% Silverman rule-of-thumb value).
h0 = 3.488946569653070;

% Numerical optimisation to find bandwidth h (using S1 values as input)
h_S1_ml = fminunc(@(h) -objfcn_S1_ml(h), h0);

%=====
% NOTE
%=====
% Currently using the average of h_S1_ml and h_S2_ml.
% Because of the random permutation applied on the wage vector at the
% beginning, MLE estimations using S1 and S2 subsamples always change. The
% reference estimates used for the following density plot are as follows:
% h_S1_ml = 0.123325277944177;
% h_S2_ml = 1.414264623487574.
%=========
% END NOTE
%=========
h = (h_S1_ml + h_S2_ml)/2;

% Creating Gaussian kernel function (input is experience)
K_exper = @(x_exper) normpdf((x_exper - exper)./h);

% Creating indicator function for female (no kernel function because gender
% is binary/discrete)
indicator_female = gender == 1;

% Creating indicator function for male (no kernel function because gender
% is binary/discrete)
indicator_male = gender == 0;

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are experience and gender (female))
m_numerator_female = @(x_exper) sum(wage.*K_exper(x_exper).*indicator_female);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are experience and gender (female))
m_denominator_female = @(x_exper) sum(K_exper(x_exper).*indicator_female);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% experience and gender(female))
m_female = @(x_exper) m_numerator_female(x_exper)./m_denominator_female(x_exper);

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are experience and gender (male))
m_numerator_male = @(x_exper) sum(wage.*K_exper(x_exper).*indicator_male);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are experience and gender (male))
m_denominator_male = @(x_exper) sum(K_exper(x_exper).*indicator_male);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% experience and gender(female))
m_male = @(x_exper) m_numerator_male(x_exper)./m_denominator_male(x_exper);

% Creating kernel conditional expectation estimation plots
fplot(m_female, [min(exper) max(exper)], 'LineWidth', 2)
hold on
fplot(m_male, [min(exper) max(exper)], 'LineWidth', 2)
grid on
xlabel('Experience');
ylabel('m(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = 'gender = female = 1';
Legend{2} = 'gender = male = 0';
legend(Legend, 'Location', 'Northeast');

hold on
title({'Kernel conditional expectation estimator', '$\hat{m}_{h}(x) \approx E[wage | experience, gender]$,', 'Using Gaussian kernel,' ,'h is the average of MLE estimates found', 'from cross-validation via splitting the sample', 'into equal-sized sub-samples S1 and S2'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\1iib_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 1iii: Use the kernel method to estimate both the conditional expectation of wage on a subset values of covariates (you pick the subset)
%=====
% NOTE
%=====
% Because we are unsure of the problem's wording, we will for now estimate
% E[wage | tenure, married].
%=========
% END NOTE
%=========
%====================================================================
% We are choosing bandwidth h based on the rule-of-thumb method. More
% information can be found at:
% https://en.wikipedia.org/wiki/Kernel_density_estimation#A_rule-of-thumb_bandwidth_estimator
% More specifically, we are using Silverman's rule-of-thumb.
%====================================================================
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Creating vector of tenure from data table
tenure = table2array(wage1(:, 5));
% Creating vector of marriage status (married = 1, not married = 0) from
% data table
married = table2array(wage1(:, 8));

% Getting sample size
n = length(wage);
% Calculating standard deviation of sample wages
sigma_hat = std(tenure);
% Calculating bandwidth h according to Silverman's rule-of-thumb and with
% respect to tenure
h = 0.9*min(sigma_hat, (iqr(tenure)/1.34))*(n^(-1/5));

% Creating Gaussian kernel function (input is tenure)
K_tenure = @(x_tenure) normpdf((x_tenure - tenure)./h);

% Creating indicator function for marriage status (no kernel function
% because gender is binary/discrete)
indicator_married = married == 1;

% Creating indicator function for marriage status (no kernel function
% because gender is binary/discrete)
indicator_notmarried = married == 0;

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are tenure and marriage status (married))
m_numerator_married = @(x_tenure) sum(wage.*K_tenure(x_tenure).*indicator_married);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are tenure and marriage status (married))
m_denominator_married = @(x_tenure) sum(K_tenure(x_tenure).*indicator_married);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% tenure and marriage status (married))
m_married = @(x_tenure) m_numerator_married(x_tenure)./m_denominator_married(x_tenure);

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are tenure and marriage status (not married))
m_numerator_notmarried = @(x_tenure) sum(wage.*K_tenure(x_tenure).*indicator_notmarried);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are tenure and marriage status (not married))
m_denominator_notmarried = @(x_tenure) sum(K_tenure(x_tenure).*indicator_notmarried);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% tenure and marriage status (not married))
m_notmarried = @(x_tenure) m_numerator_notmarried(x_tenure)./m_denominator_notmarried(x_tenure);

% Creating kernel conditional expectation estimation plots
fplot(m_married, [min(tenure) max(tenure)], 'LineWidth', 2)
hold on
fplot(m_notmarried, [min(tenure) max(tenure)], 'LineWidth', 2)
grid on
xlabel('Tenure');
ylabel('m(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = 'Marriage status = married = 1';
Legend{2} = 'Marriage status = not married = 0';
legend(Legend, 'Location', 'Southwest');

hold on
title({'Kernel conditional expectation estimator', '$\hat{m}_{h}(x) \approx E[wage | tenure, marriage status]$,', 'Using Gaussian kernel,' ,'h calculated using the Silverman rule-of-thumb.'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\1iiia_plot.png');
close(gcf);

%=======================================================
% We are choosing bandwidth h based on cross-validation.
%=======================================================
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Creating vector of tenure from data table
tenure = table2array(wage1(:, 5));
% Creating vector of tenure from data table. For cross-validation purposes,
% we are randomising the tenure vector to obtain better distributed testing
% and training sets.
tenure_rand = tenure(randperm(length(tenure)));
% Creating vector of marriage status (married = 1, not married = 0) from
% data table
married = table2array(wage1(:, 8));

% Getting sample size
n = length(tenure_rand);
% Creating subsample S1
global S1;
S1 = tenure_rand(1:(n/2), 1);
% Getting S1 size
global n1;
n1 = length(S1);
% Creating subsample S2
global S2;
S2 = tenure_rand(((n/2) + 1):n, 1);
% Getting S2 size
global n2;
n2 = length(S2);

%=====
% NOTE
%=====
% Obtaining MLE of S2 sample
%=========
% END NOTE
%=========
% Creating kernel function that takes a single S2 value as input
global K_S2_i;
K_S2_i = @(x, h) normpdf((S1 - x)./h);

% Creating density function that takes a single S2 value as input
global f_S2_i;
f_S2_i = @(x, h) sum(K_S2_i(x, h))./(n1.*h);

% Creating log-likelihood contribution of a single S2 value using above
% density function
global ll_S2_i;
ll_S2_i = @(x, h) log(f_S2_i(x, h));

% Creating objective function for MLE (takes S2 values as input for
% calculation). More information about the ll_S2 function can be found in
% ll_S2.m
objfcn_S2_ml = @(h) ll_S2(h);

% Initial bandwidth value for numerical optimisation (setting equal to
% Silverman rule-of-thumb)
h0 = 1.342884318084240;

% Numerical optimisation to find bandwidth h (using S2 values as input)
h_S2_ml = fminunc(@(h) -objfcn_S2_ml(h), h0);

%=====
% NOTE
%=====
% Obtaining MLE of S1 sample
%=========
% END NOTE
%=========
% Creating kernel function that takes a single S1 value as input
global K_S1_i;
K_S1_i = @(x, h) normpdf((S2 - x)./h);

% Creating density function that takes a single S1 value as input
global f_S1_i;
f_S1_i = @(x, h) sum(K_S1_i(x, h))./(n2.*h);

% Creating log-likelihood contribution of a single S1 value using above
% density function
global ll_S1_i;
ll_S1_i = @(x, h) log(f_S1_i(x, h));

% Creating objective function for MLE (takes S1 values as input for
% calculation). More information about the ll_S1 function can be found in
% ll_S1.m
objfcn_S1_ml = @(h) ll_S1(h);

% Initial bandwidth value for numerical optimisation (setting equal to the
% Silverman rule-of-thumb value).
h0 = 1.342884318084240;

% Numerical optimisation to find bandwidth h (using S1 values as input)
h_S1_ml = fminunc(@(h) -objfcn_S1_ml(h), h0);

%=====
% NOTE
%=====
% Currently using the average of h_S1_ml and h_S2_ml.
% Because of the random permutation applied on the wage vector at the
% beginning, MLE estimations using S1 and S2 subsamples always change. The
% reference estimates used for the following density plot are as follows:
% h_S1_ml = 0.883360242883666;
% h_S2_ml = 0.851014469418581.
%=========
% END NOTE
%=========
h = (h_S1_ml + h_S2_ml)/2;

% Creating Gaussian kernel function (input is tenure)
K_tenure = @(x_tenure) normpdf((x_tenure - tenure)./h);

% Creating indicator function for marriage status (no kernel function
% because gender is binary/discrete)
indicator_married = married == 1;

% Creating indicator function for marriage status (no kernel function
% because gender is binary/discrete)
indicator_notmarried = married == 0;

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are tenure and marriage status (married))
m_numerator_married = @(x_tenure) sum(wage.*K_tenure(x_tenure).*indicator_married);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are tenure and marriage status (married))
m_denominator_married = @(x_tenure) sum(K_tenure(x_tenure).*indicator_married);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% tenure and marriage status (married))
m_married = @(x_tenure) m_numerator_married(x_tenure)./m_denominator_married(x_tenure);

% Creating numerator of the Nadaraya–Watson kernel regression function
% (inputs are tenure and marriage status (not married))
m_numerator_notmarried = @(x_tenure) sum(wage.*K_tenure(x_tenure).*indicator_notmarried);

% Creating denominator of the Nadaraya–Watson kernel regression function
% (inputs are tenure and marriage status (not married))
m_denominator_notmarried = @(x_tenure) sum(K_tenure(x_tenure).*indicator_notmarried);

% Creating the Nadaraya–Watson kernel regression function (inputs are
% tenure and marriage status (not married))
m_notmarried = @(x_tenure) m_numerator_notmarried(x_tenure)./m_denominator_notmarried(x_tenure);

% Creating kernel conditional expectation estimation plots
fplot(m_married, [min(tenure) max(tenure)], 'LineWidth', 2)
hold on
fplot(m_notmarried, [min(tenure) max(tenure)], 'LineWidth', 2)
grid on
xlabel('Tenure');
ylabel('m(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = 'Marriage status = married = 1';
Legend{2} = 'Marriage status = not married = 0';
legend(Legend, 'Location', 'Southwest');

hold on
title({'Kernel conditional expectation estimator', '$\hat{m}_{h}(x) \approx E[wage | tenure, marriage status]$,', 'Using Gaussian kernel,' ,'h is the average of MLE estimates found', 'from cross-validation via splitting the sample', 'into equal-sized sub-samples S1 and S2'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\1iiib_plot.png');
close(gcf);

%==========================================================================