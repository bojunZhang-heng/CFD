clc; clear; close all;

%------------------------------------------------------------------------
% Load data
data_1 = load('L1_error_u_N.dat');  
data_2 = load('L1_error_u_M.dat');  
data_3 = load('L1_error_u_AS.dat');  

% Define x-coordinates
x_coor = [100, 200, 300, 400, 500]';

% Extract errors
L1_error_u_N  = data_1(:,2);
L1_error_u_M  = data_2(:,2);
L1_error_u_AS = data_3(:,2);

%------------------------------------------------------------------------
% Compute slopes
slope_N  = polyfit(log(x_coor), log(L1_error_u_N), 1);
slope_M  = polyfit(log(x_coor), log(L1_error_u_M), 1);
slope_AS = polyfit(log(x_coor), log(L1_error_u_AS), 1);

fprintf('Slope for L1_u_N: %.2f\n', slope_N(1));
fprintf('Slope for L1_u_M: %.2f\n', slope_M(1));
fprintf('Slope for L1_u_AS: %.2f\n', slope_AS(1));

%------------------------------------------------------------------------
% Create plot
figure('Position', [100, 100, 800, 600]);
hold on;

% Plot data with LaTeX-escaped legend entries
loglog(x_coor, L1_error_u_N, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, ...
       'DisplayName', '$u\_N$'); % Escaped underscore
loglog(x_coor, L1_error_u_M, 'b-s', 'LineWidth', 2, 'MarkerSize', 8, ...
       'DisplayName', '$u\_M$');
loglog(x_coor, L1_error_u_AS, 'g-^', 'LineWidth', 2, 'MarkerSize', 8, ...
       'DisplayName', '$u\_{AS}$');

% Reference line
C = L1_error_u_N(1) * x_coor(1)^2;
loglog(x_coor, C*x_coor.^(-2), 'k--', 'LineWidth', 2, ...
       'DisplayName', 'Slope = -2');

% Formatting
set(gca, 'XScale', 'log', 'YScale', 'log');
xticks(x_coor);
xlabel('Grid Resolution (N)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('L1 Error', 'Interpreter', 'latex', 'FontSize', 14);
title('Error Convergence Analysis', 'Interpreter', 'latex', 'FontSize', 16);

% Legend with LaTeX interpreter
lgd = legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'latex');
lgd.Title.String = 'Components';
grid on;

% PDF export settings
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8 6]);
saveas(gcf, 'error_convergence.pdf');
