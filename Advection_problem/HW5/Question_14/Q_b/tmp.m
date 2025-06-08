clc; clear; close all;

% load data
data_1 = readmatrix('CFL0.500.txt');  
data_2 = readmatrix('CFL1.000.txt');  
data_3 = readmatrix('CFL1.200.txt');  
%data_4 = load('L2_error_1.0.txt');  

% Extract the grid resolution N(ii)
tt = data_3(:,1);

% Extract the L1_error
L1_error1 = data_3(:,4);
L1_error2 = data_3(:,5);



% Create figure with custom size
figure('Position', [100, 100, 1400, 600]);

% The first subplot
% Custom position: [left, bottom, width, height]
%ax1 = axes('Position', [0.07, 0.15, 0.4, 0.78]); 
hold on;
set(gca, 'FontSize', 13);
plot(tt, L1_error1, 'r-o', 'LineWidth', 2, 'MarkerSize', 5);
plot(tt, L1_error2, 'b-o', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('N', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('L1\_error', 'Interpreter', 'latex', 'FontSize', 20);
title('L1\_error at $CFL = 1.5$', 'Interpreter', 'latex', 'FontSize', 20);
legend('Location', 'northeast', 'FontSize', 16, 'Interpreter', 'latex');
grid on;
hold off;

% save as pdf
set(gcf, 'PaperPositionMode', 'auto');
print('L1-error_Against_tt_3.pdf', '-dpdf', '-bestfit');

