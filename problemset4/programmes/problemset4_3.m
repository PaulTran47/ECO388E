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
% ECO388E Problem Set 4 (1, Spring 2022), 3
% Paul Le Tran, plt377
% 14 March, 2022
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
wage1 = readtable("path\to\data\wage1.csv", opts);
% Clear temporary variables
clear opts

% Returning to home directory
cd(home_dir);
%==========================================================================

%==========================================================================
%% Part 3i: Estimate the density of wage with the uniform kernel and h = \lambda*\hat{\sigma}*(n^{-1/5}), where \lambda = 0.5, 1, 1.5, 2, 2.5.
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Getting sample size
n = length(wage);
% Creating vector of \lambda values
lambda = (0.5:0.5:2.5);
% Calculating standard deviation of sample wages
sigma_hat = std(wage);
% Creating vector of bandwidths
h = lambda.*sigma_hat.*(n^(-1/5));

% Creating kernel function
K = @(x) (x - wage)./h <= 0.5 & (x - wage)./h >= -0.5;

% Creating kernel density estimator function
f = @(x) sum(K(x))./(n.*h);

% Plotting our kernel density estimator function
fplot(f, [min(wage), max(wage)], 'LineWidth', 2);
grid on
xlabel('x');
ylabel('f(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = '\lambda = 0.5';
Legend{2} = '\lambda = 1';
Legend{3} = '\lambda = 1.5';
Legend{4} = '\lambda = 2';
Legend{5} = '\lambda = 2.5';
legend(Legend, 'Location', 'best');

hold on
title({'Kernel density estimator $\hat{f}_{h}(x)$,', 'Using uniform kernel,' ,'calculated for different $\lambda$ values'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\3i_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 3ii: Estimate the density of wage with the Gaussian kernel and h = \lambda*\hat{\sigma}*(n^{-1/5}), where \lambda = 0.5, 1, 1.5, 2, 2.5.
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Getting sample size
n = length(wage);
% Creating vector of \lambda values
lambda = (0.5:0.5:2.5);
% Calculating standard deviation of sample wages
sigma_hat = std(wage);
% Creating vector of bandwidths
h = lambda.*sigma_hat.*(n^(-1/5));

% Creating kernel function
K = @(x) normpdf((x - wage)./h);

% Creating kernel density estimator function
f = @(x) sum(K(x))./(n.*h);

% Plotting our kernel density estimator function
fplot(f, [min(wage), max(wage)], 'LineWidth', 2);
grid on
xlabel('x');
ylabel('f(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = '\lambda = 0.5';
Legend{2} = '\lambda = 1';
Legend{3} = '\lambda = 1.5';
Legend{4} = '\lambda = 2';
Legend{5} = '\lambda = 2.5';
legend(Legend, 'Location', 'best');

hold on
title({'Kernel density estimator $\hat{f}_{h}(x)$,', 'Using Gaussian kernel,' ,'calculated for different $\lambda$ values'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\3ii_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 3iii: Estimate the density of wage with the Epanechnikov kernel and h = \lambda*\hat{\sigma}*(n^{-1/5}), where \lambda = 0.5, 1, 1.5, 2, 2.5.
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Getting sample size
n = length(wage);
% Creating vector of \lambda values
lambda = (0.5:0.5:2.5);
% Calculating standard deviation of sample wages
sigma_hat = std(wage);
% Creating vector of bandwidths
h = lambda.*sigma_hat.*(n^(-1/5));

% Creating indicator function component of kernel function
K_indicator = @(x) (x - wage)./h <= 1 & (x - wage)./h >= -1;

% Creating kernel function
K = @(x) 0.75*(1 - ((x - wage)./h).^2).*K_indicator(x);

% Creating kernel density estimator function
f = @(x) sum(K(x))./(n.*h);

% Plotting our kernel density estimator function
fplot(f, [min(wage), max(wage)], 'LineWidth', 2);
grid on
xlabel('x');
ylabel('f(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = '\lambda = 0.5';
Legend{2} = '\lambda = 1';
Legend{3} = '\lambda = 1.5';
Legend{4} = '\lambda = 2';
Legend{5} = '\lambda = 2.5';
legend(Legend, 'Location', 'best');

hold on
title({'Kernel density estimator $\hat{f}_{h}(x)$,', 'Using Epanechnikov kernel,' ,'calculated for different $\lambda$ values'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\3iii_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 3iv: Estimate the density of wage with the Gaussian kernel and h = 1.06*\hat{\sigma}*(n^{-r}), where r = 1/2, 1/3, 1/4, 1/5, 1/6, 1/7.
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Getting sample size
n = length(wage);
% Creating vector of r values
r = [1/2 1/3 1/4 1/5 1/6 1/7];
% Calculating standard deviation of sample wages
sigma_hat = std(wage);
% Creating vector of bandwidths
h = 1.06*sigma_hat.*(n.^(-r));

% Creating kernel function
K = @(x) normpdf((x - wage)./h);

% Creating kernel density estimator function
f = @(x) sum(K(x))./(n.*h);

% Plotting our kernel density estimator function
fplot(f, [min(wage), max(wage)], 'LineWidth', 2);
grid on
xlabel('x');
ylabel('f(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = 'r = 1/2';
Legend{2} = 'r = 1/3';
Legend{3} = 'r = 1/4';
Legend{4} = 'r = 1/5';
Legend{5} = 'r = 1/6';
Legend{6} = 'r = 1/7';
legend(Legend, 'Location', 'best');

hold on
title({'Kernel density estimator $\hat{f}_{h}(x)$,', 'Using Gaussian kernel,' ,'calculated for different $r$ values'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\3iv_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 3v: Estimate the density of wage with the Gaussian kernel and h is chosen by the cross validation method.
%=====
% NOTE
%=====
% We are choosing to split the sample rather than do LOOCV. We will split
% into two sub-samples of equal size. Specifically, we will create the
% kernel density estimator whose inputs are values in subsample S2. When
% then perform MLE to find the optimal bandwidth h.
%=========
% END NOTE
%=========
% Creating vector of wages from data table. For cross-validation purposes,
% we are randomising the wage vector to obtain better distributed testing
% and training sets.
wage = table2array(wage1(:, 2));
wage_rand = wage(randperm(length(wage)));
% Getting sample size
n = length(wage_rand);
% Creating subsample S1
global S1;
S1 = wage_rand(1:263, 1);
% Getting S1 size
global n1;
n1 = length(S1);
% Creating subsample S2
global S2;
S2 = wage_rand(264:n, 1);
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

% Initial bandwidth value for numerical optimisation
h0 = 0.5;

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

% Initial bandwidth value for numerical optimisation
h0 = 0.5;

% Numerical optimisation to find bandwidth h (using S1 values as input)
h_S1_ml = fminunc(@(h) -objfcn_S1_ml(h), h0);

% Setting our bandwidth to be the MLE result
%=====
% NOTE
%=====
% Currently using the average of h_S1_ml and h_S2_ml.
% Because of the random permutation applied on the wage vector at the
% beginning, MLE estimations using S1 and S2 subsamples always change. The
% reference estimates used for the following density plot are as follows:
% h_S1_ml = 0.6022;
% h_S2_ml = 0.5227.
%=========
% END NOTE
%=========
h = (h_S1_ml + h_S2_ml)/2;

% Creating general kernel function
K = @(x) normpdf((x - wage)./h);

% Creating kernel density estimator function
f = @(x) sum(K(x))./(n.*h);

% Plotting our kernel density estimator function
fplot(f, [min(wage), max(wage)], 'LineWidth', 2);
grid on
xlabel('x');
ylabel('f(x)');
hold on

hold on
title({'Kernel density estimator $\hat{f}_{h}(x)$,', 'Using Gaussian kernel,' ,'h is the average of MLE estimates found', 'from cross-validation via splitting the sample', 'into equal-sized sub-samples S1 and S2'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\3v_plot.png');
close(gcf);
%==========================================================================