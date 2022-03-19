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
% ECO388E Problem Set 4 (1, Spring 2022), 2
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
%% Part 2i: Plot the empirical CDF function of wage estimated with randomly picking 10 observations from the sample
% Picking 10 random observations
obs10 = table2array(datasample(wage1(:, 2), 10));

% Plotting empirical CDF of wage
set(cdfplot(obs10), 'LineWidth', 2);
hold on
grid on
xlabel('x');
ylabel('F(x)');
title({'Empirical CDF of wage, calculated with', '10 random observations'});

saveas(gcf, 'path\to\graphics\2i_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 2ii: Plot the empirical CDF function of wage estimated with randomly picking 100 observations from the sample
% Picking 100 random observations
obs100 = table2array(datasample(wage1(:, 2), 100));

% Plotting the empirical CDFs of wage calculated using 10 and 100 random
% observations from sample
set(cdfplot(obs10), 'LineWidth', 2);
hold on
set(cdfplot(obs100), 'LineWidth', 2);
hold on
grid on
xlabel('x');
ylabel('F(x)');
legend('10 random observations', '100 random observations', 'Location', 'best');
hold on
title({'Empirical CDFs of wage, calculated with', '10 and 100 random observations respectively'});

saveas(gcf, 'path\to\graphics\2ii_plot.png');
close(gcf);
%==========================================================================
%% Part 2iii: Plot the empirical CDF function of wage estimated with the full sample
% Creating full sample of wages into array
obs = table2array(wage1(:, 2));

% Plotting the empirical CDFs of wage calculated using 10 and 100 random
% samples from the sample, and the full sample
set(cdfplot(obs10), 'LineWidth', 2);
hold on
set(cdfplot(obs100), 'LineWidth', 2);
hold on
set(cdfplot(obs), 'LineWidth', 2);
hold on
grid on
xlabel('x');
ylabel('F(x)');
legend('10 random observations', '100 random observations', 'Full sample', 'Location', 'best');
hold on
title({'Empirical CDFs of wage, calculated with', '10 and 100 random observations, and full sample'});

saveas(gcf, 'path\to\graphics\2iii_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 2iv: Plot the empirical CDF function of wage estimated using the estimator, \widetilde{F}(x)
%=====
% NOTE
%=====
% \widetilde{F}(x) = \frac{1}{n} \sum_{i = 1}^{n} \Phi(\frac{x - X_{i}}{h}) 
%=========
% END NOTE
%=========

% Getting number of all observations
n = length(obs);

% Creating a vector of bandwidth estimates
h = [n^(-1/5) n^(-1/4) n^(-1/3) n^(-1/2)];

% Creating our estimate function for different bandwidths
F_tilde = @(x) sum(normcdf((x - obs)./h))./n;

% Plotting the empirical CDFs of wage calculated using 10 and 100 random
% samples from the sample, the full sample, and the estimator
% \wildtilde{F}(x) for different bandwidths
set(cdfplot(obs10), 'LineWidth', 2);
hold on
set(cdfplot(obs100), 'LineWidth', 2);
hold on
set(cdfplot(obs), 'LineWidth', 2);
hold on
fplot(F_tilde, [0, 25], 'LineWidth', 2);
hold on
grid on
xlabel('x');
ylabel('F(x)');
legend('10 random observations', '100 random observations', 'Full sample', 'h = 100^{(-1/5)}', 'h = 100^{(-1/4)}', 'h = 100^{(-1/3)}', 'h = 100^{(-1/2)}', 'Location', 'best');
hold on
title({'Estimate function $\widetilde{F}(x)$'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\2iv_plot.png');
close(gcf);
%==========================================================================