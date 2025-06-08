clc; clear all;

% Load the data from the Fortran output file

data = load("HW1_1_result.txt");

% Extract M values (first column) and P values (second column)
M = data(:,1);
P = data(:,2);

% Ploi the data
figure;
plot(M, P, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 5);
xlabel('M (Number of Heads)');
ylabel('P (Probability)');
title('Binomial Probability Distribution');
grid on;

% Save the plot as an image (optional)
saveas(gcf, 'binomial_plot.png');
