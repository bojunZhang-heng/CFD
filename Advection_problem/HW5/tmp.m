clc; clear; close all;

% load data
data_1 = readmatrix('t1.txt', 'NumHeaderLines', 1);  
data_2 = readmatrix('t2.txt', 'NumHeaderLines', 1);  
data_3 = readmatrix('t3.txt', 'NumHeaderLines', 1);  
data_4 = readmatrix('t4.txt', 'NumHeaderLines', 1);  

% Extract the grid resolution N(ii)
x_coor = data_4(:, 2);

% Extract the numberical solutiuon
u_N1 = data_4(:, 3);
u_N2 = data_4(:, 4);
u_N3 = data_4(:, 5);
u    = data_4(:, 6);



% Create figure with custom size
figure('Position', [100, 100, 1400, 600]);

% The first subplot
% Custom position: [left, bottom, width, height]
%ax1 = axes('Position', [0.07, 0.15, 0.4, 0.78]); 
hold on;
set(gca, 'FontSize', 13);
plot(x_coor, u_N1, 'c-^', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Lax-Wendroff');
plot(x_coor, u_N2, 'm-.s', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'van Leer limiter');
plot(x_coor, u_N3, 'b:x', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'SUPERBEE limiter');
plot(x_coor, u,    'g*'  , 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Analytical solu');
xlabel('x\_coor', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('solution', 'Interpreter', 'latex', 'FontSize', 20);
title('at $tt = 2.0$', 'Interpreter', 'latex', 'FontSize', 20);
legend('Location', 'northeast', 'FontSize', 8, 'Interpreter', 'latex');
grid on;
hold off;

% save as pdf
set(gcf, 'PaperPositionMode', 'auto');
print('t4.pdf', '-dpdf', '-bestfit');

