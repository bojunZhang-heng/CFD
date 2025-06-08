clc; clear; close all;

% load data
data_1 = readmatrix('t1.txt', 'NumHeaderLines', 1);  
data_2 = readmatrix('t2.txt', 'NumHeaderLines', 1);  
data_3 = readmatrix('t3.txt', 'NumHeaderLines', 1);  
data_4 = readmatrix('TV.txt', 'NumHeaderLines', 1);  

% Extract the grid resolution N(ii)
tt = data_4(:, 1);

% Extract the numberical solutiuon
TV1 = data_4(:, 2);
TV2 = data_4(:, 3);
TV3 = data_4(:, 4);



% Create figure with custom size
figure('Position', [100, 100, 1400, 600]);

% The first subplot
% Custom position: [left, bottom, width, height]
%ax1 = axes('Position', [0.07, 0.15, 0.4, 0.78]); 
hold on;
set(gca, 'FontSize', 13);
plot(tt, TV1, 'c-^', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Lax-Wendroff');
plot(tt, TV2, 'm-.s', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'van Leer limiter');
plot(tt, TV3, 'b:x', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'SUPERBEE limiter');
xlabel('x\_coor', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('solution', 'Interpreter', 'latex', 'FontSize', 20);
title('$TV$', 'Interpreter', 'latex', 'FontSize', 20);
legend('Location', 'northeast', 'FontSize', 14, 'Interpreter', 'latex');
grid on;
hold off;

% save as pdf
set(gcf, 'PaperPositionMode', 'auto');
print('TV.pdf', '-dpdf', '-bestfit');

