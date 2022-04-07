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
% ECO388E Problem Set 5 (2, Spring 2022), 2
% Paul Le Tran, plt377
% 5 April, 2022
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
%% Part 2i: Use the K-nearest neighbour approach to estimate conditional expectation of wage on education and gender
%=====
% NOTE
%=====
% Observe that our X variables are education and gender (latter is binary),
% the K-nearest neighbour approach will revolve around these two variables.
% More specifically, what we will do is find the K-nearest neighbours for
% female and male.

% We will choose K = \lambda*\sqrt(n), where we pick \lambda by CV.
%=========
% END NOTE
%=========
%================================
% Performing CV to obtain \lambda
%================================
% Getting total sample size
n = height(wage1);
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Creating vector of education from data table
edu = table2array(wage1(:, 3));
% Creating vector of gender (female = 1, female = 0) from data table
gender = table2array(wage1(:, 7));

% Putting education, gender, and wage into a single matrix for shuffling
% purposes
edu_gender_wage = [edu gender wage];
edu_gender_wage_rand = edu_gender_wage(randperm(size(edu_gender_wage, 1)), :);
% Creating training and testing subsets by a 80-20 split rule.
global train test;
train = edu_gender_wage_rand(1:round(0.8*length(edu_gender_wage_rand)), :);
test = edu_gender_wage_rand(round(0.8*length(edu_gender_wage_rand)):end, :);

% Calculating the true mean wage conditional on education and gender
% of the test subset.
% Initialising matrix that will house all unique education values and
% associatied mean wage from the test subset.
global test_unique;
test_unique = unique(test(:, 1));
test_true = zeros(length(test_unique), 2);
test_true(:, 1) = test_unique;

for i = 1:length(test_unique)
  disp(i);
  test_true(i, 2) = mean(test(test(:,1) == test_true(i, 1), 2));
end

% Creating K (function of lambda)
global K;
K = @(lambda) lambda.*sqrt(n);

% Creating objective function. More information about the specifics done to
% calculate the MSE can be found in mse.m
obj_func = @(lambda) mse(lambda);

% Initial \lambda value for numerical optimisation. We choose the starting
% value based on how similar the resulting plots look to problem 1 results
% as a test.
lambda0 = 2.2;

% Numerical optimisation to find bandwidth \lambda
%lambda_star = fminsearch(@(lambda) obj_func(lambda), lambda0);
%=======
% ANSWER
%=======
% The optimal \lambda^{*} is dependent on the initial dataset split. For
% reproducibility purposes, we will save \lambda^{*} = 2.4200 manually and
% comment out the numerical optimisation step.
%===========
% END ANSWER
%===========
lambda_star = 2.42;

%================================
% Performing KNN for only females
%================================
% Filtering table to only display females
wage1_female = wage1(wage1.female == 1,:);
% Creating vector of wages from data table
wage = table2array(wage1_female(:, 2));
% Creating vector of education from data table
edu = table2array(wage1_female(:, 3));

% Getting total sample size
n = height(wage1);
% Defining K, which is the number of neighbours we are considering.
K_star = round(lambda_star.*sqrt(n));

% Declaring matrix that will house estimate of E[wage | education, female]
m_knn_female_edu = zeros(max(wage1.educ) + 1, 1);

% Performing KNN search to estiamte E[wage | education, female]
% Declaring vector that will hold associated estimate of
% E[wage | education, female].
m_knn_female = zeros(max(wage1.educ) + 1, 1); 
for j = 1:max(wage1.educ) + 1
  disp(j - 1);
  Idx = knnsearch(edu, j - 1, 'K', K_star);
  
  % Finding average of corresponding wage values to these neighbours.
  m_knn_female(j) = mean(wage(Idx, :));
end
  
% Adding column of average corresponding wages to main matrix
m_knn_female_edu(:, 2) = m_knn_female;

%==============================
% Performing KNN for only males
%==============================
% Filtering table to only display males
wage1_male = wage1(wage1.female == 0,:);
% Creating vector of wages from data table
wage = table2array(wage1_male(:, 2));
% Creating vector of education from data table
edu = table2array(wage1_male(:, 3));

% Getting total sample size
n = height(wage1);
% Defining K, which is the number of neighbours we are considering.
K_star = round(lambda_star.*sqrt(n));

% Declaring matrix that will house estimate of E[wage | education, male]
m_knn_male_edu = zeros(max(wage1.educ) + 1, 1);

% Performing KNN search to estiamte E[wage | education, male]
% Declaring vector that will hold associated estimate of
% E[wage | education, male].
m_knn_male = zeros(max(wage1.educ) + 1, 1); 
for j = 1:max(wage1.educ) + 1
  disp(j - 1);
  Idx = knnsearch(edu, j - 1, 'K', K_star);
  
  % Finding average of corresponding wage values to these neighbours.
  m_knn_male(j) = mean(wage(Idx, :));
end
  
% Adding column of average corresponding wages to main matrix
m_knn_male_edu(:, 2) = m_knn_male;

%=================
% Plotting results
%=================
% Creating KNN-search conditional expectation estimation plots for males
plot((0:1:max(wage1.educ))', m_knn_female_edu(:, 2), 'LineWidth', 2);
hold on;
plot((0:1:max(wage1.educ))', m_knn_male_edu(:, 2), 'LineWidth', 2);
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
title({'KNN-search conditional expectation estimator', '$\hat{m}_{h}(x) \approx E[wage | education, gender]$,', 'calculated K (and $\lambda$) using cross-validation'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\2i_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 2ii: Use the K-nearest neighbour approach to estimate conditional expectation of wage on experience and gender
%=====
% NOTE
%=====
% Observe that our X variables are experience and gender (latter is
% binary), the K-nearest neighbour approach will revolve around these two
% variables. More specifically, what we will do is find the K-nearest
% neighbours for female and male.

% We will choose K = \lambda*\sqrt(n), where we pick \lambda by CV.
%=========
% END NOTE
%=========
%================================
% Performing CV to obtain \lambda
%================================
% Getting total sample size
n = height(wage1);
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Creating vector of experience from data table
exper = table2array(wage1(:, 4));
% Creating vector of gender (female = 1, female = 0) from data table
gender = table2array(wage1(:, 7));

% Putting experience, gender, and wage into a single matrix for shuffling
% purposes
exper_gender_wage = [exper gender wage];
exper_gender_wage_rand = exper_gender_wage(randperm(size(exper_gender_wage, 1)), :);
% Creating training and testing subsets by a 80-20 split rule.
global train test;
train = exper_gender_wage_rand(1:round(0.8*length(exper_gender_wage_rand)), :);
test = exper_gender_wage_rand(round(0.8*length(exper_gender_wage_rand)):end, :);

% Calculating the true mean wage conditional on experience and gender
% of the test subset.
% Initialising matrix that will house all unique experience values and
% associatied mean wage from the test subset.
global test_unique;
test_unique = unique(test(:, 1));
test_true = zeros(length(test_unique), 2);
test_true(:, 1) = test_unique;

for i = 1:length(test_unique)
  disp(i);
  test_true(i, 2) = mean(test(test(:,1) == test_true(i, 1), 2));
end

% Creating K (function of lambda)
global K;
K = @(lambda) lambda.*sqrt(n);

% Creating objective function. More information about the specifics done to
% calculate the MSE can be found in mse.m
obj_func = @(lambda) mse(lambda);

% Initial \lambda value for numerical optimisation. We choose the starting
% value based on how similar the resulting plots look to problem 1 results
% as a test.
lambda0 = 2.2;

% Numerical optimisation to find bandwidth \lambda
%lambda_star = fminsearch(@(lambda) obj_func(lambda), lambda0);
%=======
% ANSWER
%=======
% The optimal \lambda^{*} is dependent on the initial dataset split. For
% reproducibility purposes, we will save \lambda^{*} = 2.365 manually and
% comment out the numerical optimisation step.
%===========
% END ANSWER
%===========
lambda_star = 2.365;

%================================
% Performing KNN for only females
%================================
% Filtering table to only display females
wage1_female = wage1(wage1.female == 1,:);
% Creating vector of wages from data table
wage = table2array(wage1_female(:, 2));
% Creating vector of experience from data table
exper = table2array(wage1_female(:, 4));

% Getting total sample size
n = height(wage1);
% Defining K, which is the number of neighbours we are considering.
K_star = round(lambda_star.*sqrt(n));

% Declaring matrix that will house estimate of E[wage | experience, female]
m_knn_female_exper = zeros(max(wage1.exper) + 1, 1);

% Performing KNN search to estiamte E[wage | experience, female]
% Declaring vector that will hold associated estimate of
% E[wage | experience, female].
m_knn_female = zeros(max(wage1.exper) + 1, 1); 
for j = 1:max(wage1.exper) + 1
  disp(j - 1);
  Idx = knnsearch(exper, j - 1, 'K', K_star);
  
  % Finding average of corresponding wage values to these neighbours.
  m_knn_female(j) = mean(wage(Idx, :));
end
  
% Adding column of average corresponding wages to main matrix
m_knn_female_exper(:, 2) = m_knn_female;

%==============================
% Performing KNN for only males
%==============================
% Filtering table to only display males
wage1_male = wage1(wage1.female == 0,:);
% Creating vector of wages from data table
wage = table2array(wage1_male(:, 2));
% Creating vector of experience from data table
exper = table2array(wage1_male(:, 4));

% Getting total sample size
n = height(wage1);
% Defining K, which is the number of neighbours we are considering.
K_star = round(lambda_star.*sqrt(n));

% Declaring matrix that will house estimate of E[wage | experience, male]
m_knn_male_exper = zeros(max(wage1.exper) + 1, 1);

% Performing KNN search to estiamte E[wage | experience, male]
% Declaring vector that will hold associated estimate of
% E[wage | experience, male].
m_knn_male = zeros(max(wage1.exper) + 1, 1); 
for j = 1:max(wage1.exper) + 1
  disp(j - 1);
  Idx = knnsearch(exper, j - 1, 'K', K_star);
  
  % Finding average of corresponding wage values to these neighbours.
  m_knn_male(j) = mean(wage(Idx, :));
end
  
% Adding column of average corresponding wages to main matrix
m_knn_male_exper(:, 2) = m_knn_male;

%=================
% Plotting results
%=================
% Creating KNN-search conditional expectation estimation plots for males
plot((0:1:max(wage1.exper))', m_knn_female_exper(:, 2), 'LineWidth', 2);
hold on;
plot((0:1:max(wage1.exper))', m_knn_male_exper(:, 2), 'LineWidth', 2);
grid on
xlabel('Experience');
ylabel('m(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = 'gender = female = 1';
Legend{2} = 'gender = male = 0';
legend(Legend, 'Location', 'best');

hold on
title({'KNN-search conditional expectation estimator', '$\hat{m}_{h}(x) \approx E[wage | experience, gender]$,', 'calculated K (and $\lambda$) using cross-validation'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\2ii_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 2iii: Use the K-nearest neighbour approach to estimate conditional expectation of wage on a subset values of covariates (you pick the subset)
%=====
% NOTE
%=====
% Observe that our X variables are tenure and marriage status (latter is
% binary), the K-nearest neighbour approach will revolve around these two
% variables. More specifically, what we will do is find the K-nearest
% neighbours for married and not married.

% We will choose K = \lambda*\sqrt(n), where we pick \lambda by CV.
%=========
% END NOTE
%=========
%================================
% Performing CV to obtain \lambda
%================================
% Getting total sample size
n = height(wage1);
% Creating vector of wages from data table
wage = table2array(wage1(:, 2));
% Creating vector of tenure from data table
tenure = table2array(wage1(:, 5));
% Creating vector of marriage status (married = 1, not married = 0) from
% data table
married = table2array(wage1(:, 8));

% Putting tenure, marriage status, and wage into a single matrix for
% shuffling purposes
tenure_married_wage = [tenure married wage];
tenure_married_wage_rand = tenure_married_wage(randperm(size(tenure_married_wage, 1)), :);
% Creating training and testing subsets by a 80-20 split rule.
global train test;
train = tenure_married_wage_rand(1:round(0.8*length(tenure_married_wage_rand)), :);
test = tenure_married_wage_rand(round(0.8*length(tenure_married_wage_rand)):end, :);

% Calculating the true mean wage conditional on tenure and marriage status
% of the test subset.
% Initialising matrix that will house all unique tenure values and
% associatied mean wage from the test subset.
global test_unique;
test_unique = unique(test(:, 1));
test_true = zeros(length(test_unique), 2);
test_true(:, 1) = test_unique;

for i = 1:length(test_unique)
  disp(i);
  test_true(i, 2) = mean(test(test(:,1) == test_true(i, 1), 2));
end

% Creating K (function of lambda)
global K;
K = @(lambda) lambda.*sqrt(n);

% Creating objective function. More information about the specifics done to
% calculate the MSE can be found in mse.m
obj_func = @(lambda) mse(lambda);

% Initial \lambda value for numerical optimisation. We choose the starting
% value based on how similar the resulting plots look to problem 1 results
% as a test.
lambda0 = 0.47;

% Numerical optimisation to find bandwidth \lambda
%lambda_star = fminsearch(@(lambda) obj_func(lambda), lambda0);
%=======
% ANSWER
%=======
% The optimal \lambda^{*} is dependent on the initial dataset split. For
% reproducibility purposes, we will save \lambda^{*} = 0.4465 manually and
% comment out the numerical optimisation step.
%===========
% END ANSWER
%===========
lambda_star = 0.4465;

%================================
% Performing KNN for only married
%================================
% Filtering table to only display married
wage1_married = wage1(wage1.married == 1,:);
% Creating vector of wages from data table
wage = table2array(wage1_married(:, 2));
% Creating vector of tenure from data table
tenure = table2array(wage1_married(:, 5));

% Getting total sample size
n = height(wage1);
% Defining K, which is the number of neighbours we are considering.
K_star = round(lambda_star.*sqrt(n));

% Declaring matrix that will house estimate of E[wage | tenure, married]
m_knn_married_tenure = zeros(max(wage1.tenure) + 1, 1);

% Performing KNN search to estiamte E[wage | tenure, married]
% Declaring vector that will hold associated estimate of
% E[wage | tenure, married].
m_knn_married = zeros(max(wage1.tenure) + 1, 1); 
for j = 1:max(wage1.tenure) + 1
  disp(j - 1);
  Idx = knnsearch(tenure, j - 1, 'K', K_star);
  
  % Finding average of corresponding wage values to these neighbours.
  m_knn_married(j) = mean(wage(Idx, :));
end
  
% Adding column of average corresponding wages to main matrix
m_knn_married_tenure(:, 2) = m_knn_married;

%================================
% Performing KNN for only not married
%================================
% Filtering table to only display not married
wage1_notmarried = wage1(wage1.married == 0,:);
% Creating vector of wages from data table
wage = table2array(wage1_notmarried(:, 2));
% Creating vector of tenure from data table
tenure = table2array(wage1_notmarried(:, 5));

% Getting total sample size
n = height(wage1);
% Defining K, which is the number of neighbours we are considering.
K_star = round(lambda_star.*sqrt(n));

% Declaring matrix that will house estimate of E[wage | tenure, notmarried]
m_knn_notmarried_tenure = zeros(max(wage1.tenure) + 1, 1);

% Performing KNN search to estiamte E[wage | tenure, notmarried]
% Declaring vector that will hold associated estimate of
% E[wage | tenure, notmarried].
m_knn_notmarried = zeros(max(wage1.tenure) + 1, 1); 
for j = 1:max(wage1.tenure) + 1
  disp(j - 1);
  Idx = knnsearch(tenure, j - 1, 'K', K_star);
  
  % Finding average of corresponding wage values to these neighbours.
  m_knn_notmarried(j) = mean(wage(Idx, :));
end
  
% Adding column of average corresponding wages to main matrix
m_knn_notmarried_tenure(:, 2) = m_knn_notmarried;

%=================
% Plotting results
%=================
% Creating KNN-search conditional expectation estimation plots
plot((0:1:max(wage1.tenure))', m_knn_married_tenure(:, 2), 'LineWidth', 2);
hold on;
plot((0:1:max(wage1.tenure))', m_knn_notmarried_tenure(:, 2), 'LineWidth', 2);
grid on
xlabel('Tenure');
ylabel('m(x)');
hold on

% Creating legend
Legend = cell(2,1);
Legend{1} = 'Married = 1';
Legend{2} = 'Not married = 0';
legend(Legend, 'Location', 'best');

hold on
title({'KNN-search conditional expectation estimator', '$\hat{m}_{h}(x) \approx E[wage | tenure, marriage status]$,', 'calculated K (and $\lambda$) using cross-validation'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\2iii_plot.png');
close(gcf);
%==========================================================================