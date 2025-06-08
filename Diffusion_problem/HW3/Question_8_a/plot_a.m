    clc; clear; close all;

%------------------------------------------------------------------------
% Load data
data_1 = load('L1_error_0.2.txt');  
data_2 = load('L1_error_1.0.txt');  
data_3 = load('L2_error_0.2.txt');  
data_4 = load('L2_error_1.0.txt');  

% Extract grid resolution 
N = data_1(:,1);

% Extract error
L1_error_1 = data_1(:,3);
L1_error_2 = data_2(:,3);
L2_error_1 = data_3(:,3);
L2_error_2 = data_4(:,3);

%------------------------------------------------------------------------
% Compute the slopes(order of accuracy)
slope_L1_1 = polyfit(log(N), log(L1_error_1), 1);
slope_L1_2 = polyfit(log(N), log(L1_error_2), 1);
slope_L2_1 = polyfit(log(N), log(L2_error_1), 1);
slope_L2_2 = polyfit(log(N), log(L2_error_2), 1);

fprintf('Order of accuracy (slope) for L1 (τ=0.2): %.2f\n', slope_L1_1(1));
fprintf('Order of accuracy (slope) for L1 (τ=1.0): %.2f\n', slope_L1_2(1));
fprintf('Order of accuracy (slope) for L2 (τ=0.2): %.2f\n', slope_L2_1(1));
fprintf('Order of accuracy (slope) for L2 (τ=1.0): %.2f\n', slope_L2_2(1));

%------------------------------------------------------------------------
% Plot in separate subplots
figure('Position', [100, 100, 1000, 800]); % Adjust figure size

% First plot
subplot(1,2,1);
hold on;
loglog(N, L1_error_1, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'L1');
loglog(N, L2_error_1, 'b-s', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'L2');
xlabel('grid resolution(N)'); ylabel('error');
title('$$\tau = 0.2 \quad$$', ...
      'Interpreter', 'latex', 'FontSize', 12);
legend('show', 'Location', 'best');

grid on;
% Change x-coordinates
xticks([8 16 32 64 128]); 
set(gca, 'XScale', 'log', 'YScale', 'log');  % Ensure log-log scale
hold off;

% Second plot
subplot(1,2,2);
hold on;
loglog(N, L1_error_2, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'L1');
loglog(N, L2_error_2, 'b-s', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'L2');
xlabel('grid resolution(N)'); ylabel('error');
title('$$\tau = 1.0 \quad$$', ...
      'Interpreter', 'latex', 'FontSize', 12);
legend('show', 'Location', 'best');

grid on;
% Change x-coordinates
xticks([8 16 32 64 128]); 
set(gca, 'XScale', 'log', 'YScale', 'log');  % Ensure log-log scale
hold off;

% Save as PDF
saveas(gcf, 'error_profile.pdf');  







