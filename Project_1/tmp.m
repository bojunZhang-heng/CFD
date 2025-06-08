% Plot solution
clc; clear; close all;

% Load data
data_1 = load('solution.txt');  
x_coor = data_1(:,1);
u      = data_1(:,2);
u_M    = data_1(:,3);
u_N    = data_1(:,4);
%u_AS   = data_1(:,5);

% Create figure
figure('Position', [100, 100, 1200, 500]);
ax1 = axes('Position', [0.08, 0.15, 0.85, 0.75]); 
hold on;

% Plot with LaTeX-safe names
plot(x_coor, u  , 'r-v', 'LineWidth', 2, 'MarkerSize', 3, 'DisplayName', '$u$');
plot(x_coor, u_M, 'k--p', 'LineWidth', 2, 'MarkerSize', 3, 'DisplayName', '$u_M$');
plot(x_coor, u_N, 'b-.d', 'LineWidth', 2, 'MarkerSize', 3, 'DisplayName', '$u_N$');
%plot(x_coor, u_AS, 'm-.', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', '$u_N$');

% Labels & Title
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$u(x,t)$', 'Interpreter', 'latex', 'FontSize', 24);
title('Solution at $t = 0.5$', 'Interpreter', 'latex', 'FontSize', 18);

% Legend settings 
lgd = legend('Location', 'best', 'FontSize', 10, 'Interpreter', 'latex');
lgd.Title.String = 'Legend Title';
lgd.Title.Interpreter = 'latex';
grid on;

% Fix PDF size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 12 5]); % 12" wide x 5" tall
saveas(gcf, 'solution_profile_N20.pdf');
