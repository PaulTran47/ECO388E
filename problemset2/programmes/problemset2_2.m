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

% ECO388E Problem Set 2, 2
% Paul Le Tran, plt377
% 20 October, 2021
%==========================================================================

%==========================================================================
%% Model info
% y_i = I(theta1 + theta2*x1_i + theta3*x2_i + e_i > 0);
% x2_i = theta4 + theta5*x1_i + theta6*z_i + eta_i
%==========================================================================

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'D:\Documents\PhD\2021-2022\ECO388E\problem_sets\problemset2\programmes';
data_dir = 'D:\Documents\PhD\2021-2022\ECO388E\problem_sets\problemset2\data';

cd(home_dir);

%% Loading in data
cd(data_dir);
data1 = importdata('ps2.dat');

% Declaring variables
y_i = data1(:, 1);
x1_i = data1(:, 2);
x2_i = data1(:, 3);
z_i = data1(:, 4);

% Getting sample size
% Variable is global so it can be accessed in the separate workspaces of
% called functions
global N;
N = length(y_i);

clear data1;
cd(home_dir);
%==========================================================================

%==========================================================================
%% Part 2d: Binary logit model (MLE)
%==========================================================================

%==========================================================================
%% Questions asked in problemset 2_2 that don't involve code
% Part 2a:
% Recall that x2_i is women's education level (in years). One economic
% reason how x2_i could be correlated with e_i is that e_i could be
% capturing the wage a woman wants from working. This is because women who
% want a higher wage may choose more education. Therefore, if women make
% this decision simultaneously regarding the two variables, you would
% expect the correlation between x2_i and e_i to put upward bias on
% coefficient theta3.

% Part 2b: NO QUESTION TO ANSWER. THIS PART IS SETUP FOR PART 2c.

% Part 2c: 
% From part 2b, we assume that a woman's educational level depends on x1_i,
% z_i, and e_i. From the system of two equations, we assume that
% unobservables e_i and eta_i are jointly normal and are independent of
% x1_i and z_i. Observe that cov(e_i, eta_i) = rho = E[eta_i*e_i] -
% E[eta_i]*E[e_i] = E[eta_i*e_i], by the jointly normal assumption.

% Rho relates to the possible endogeneity of x2_i through expectations.
% Specifically, with the given system of two equations, we know that
% E[x2_i*e_i] = E[eta_i*e_i] = rho. Therefore, we can see that if x2_i is
% exogenous, we must have rho = 0. But since x2_i is assumed to be
% endogenous, rho != 0. Additionally, I expect rho > 0. This is because the
% educational level of one's parents (z_i) should be positively correlated
% with x2_i. This is because a parent often serves as a role model.

% Finally, the assumption that z_i does not directly determine y_i
% (i.e., z_i affects y_i only through x2_i) is exclusion restriction. This
% is one of the three principles that makes for a "good"/"robust"
% instrument in IV estimation. Without the exclusion restriction, theta3
% would not be identified non-parametrically.
%==========================================================================