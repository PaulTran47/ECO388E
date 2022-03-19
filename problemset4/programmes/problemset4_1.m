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
% ECO388E Problem Set 4 (1, Spring 2022), 1
% Paul Le Tran, plt377
% 14 March, 2022
%==========================================================================

%==========================================================================
%% Setting up workspace
clear all;
close all;
clc;

home_dir = 'path\to\programmes';

cd(home_dir);
%==========================================================================

%==========================================================================
%% Part 1: Draw 100 i.i.d. observations from the standard normal distribution
% These observations will be used for all parts of problem 1.
n = 100;
obs = normrnd(0, 1, [n, 1]);
%==========================================================================

%==========================================================================
%% Part 1i: Plot the CDF function.
x = linspace(min(obs), max(obs));
plot(x, evcdf(x, 0, 1), 'LineWidth', 2);
hold on
grid on
xlabel('x');
ylabel('F(x)');
hold on
title({'Theoretical standard normal CDF'});

saveas(gcf, 'path\to\graphics\1i_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 1ii: Plot the empirical CDF function
% We choose to plot the theoretical and empirical CDF functions together
% for easier comparison
set(cdfplot(obs), 'LineWidth', 2);
hold on
x = linspace(min(obs), max(obs));
plot(x, evcdf(x, 0, 1), 'LineWidth', 2);
hold on
grid on
xlabel('x');
ylabel('F(x)');
legend('Empirical CDF', 'Theoretical CDF', 'Location', 'best');
hold on
title({'Theoretical and empirical CDFs'});

saveas(gcf, 'path\to\graphics\1ii_plot.png');
close(gcf);
%==========================================================================

%==========================================================================
%% Part 1iii: Plot the provided estimate as a function of x, \widetilde{F}(x)
%=====
% NOTE
%=====
% \widetilde{F}(x) = \frac{1}{n} \sum_{i = 1}^{n} \Phi(\frac{x - X_{i}}{h}) 
%=========
% END NOTE
%=========

% Creating a vector of bandwidth estimates
h = [n^(-1/5) n^(-1/4) n^(-1/3) n^(-1/2)];
 
% Creating our estimate function for different bandwidths
F_tilde = @(x) sum(normcdf((x - obs)./h))./n;

% Plotting estimate function for different bandwidths
fplot(F_tilde, [-4, 4], 'LineWidth', 2);
grid on
xlabel('x');
ylabel('F(x)');
legend('h = 100^{(-1/5)}', 'h = 100^{(-1/4)}', 'h = 100^{(-1/3)}', 'h = 100^{(-1/2)}', 'Location', 'best');
hold on
title({'Estimate function $\widetilde{F}(x)$'}, 'Interpreter','LaTeX');

saveas(gcf, 'path\to\graphics\1iii_plot.png');
close(gcf);
%==========================================================================