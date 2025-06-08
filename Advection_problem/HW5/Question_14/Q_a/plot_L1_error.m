clc; clear; close all;

% load data
data_1 = load('fort.70');  
%data_2 = load('t2.txt');  
%data_3 = load('L2_error_0.2.txt');  
%data_4 = load('L2_error_1.0.txt');  

% Extract the grid resolution N(ii)
N = data_1(:,1);

% Extract the L1_error
L1_error = data_1(:,2);



% Create figure with custom size
figure('Position', [100, 100, 1400, 600]);

% The first subplot
% Custom position: [left, bottom, width, height]
%ax1 = axes('Position', [0.07, 0.15, 0.4, 0.78]); 
hold on;
plot(N, L1_error, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Numerical');
set(gca, 'FontSize', 13);
xlabel('N', 'Interpreter', 'latex', 'FontSize', 28);
ylabel('L1\_error', 'Interpreter', 'latex', 'FontSize', 28);
title('L1\_error at $t = 1.5$', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'northeast', 'FontSize', 20, 'Interpreter', 'latex');
grid on;
hold off;

% save as pdf
set(gcf, 'PaperPositionMode', 'auto');
print('L1-error_Against_N.pdf', '-dpdf', '-bestfit');

