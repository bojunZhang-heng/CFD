clc; clear; close all;

% load data
data_1 = load('t1.txt');  
data_2 = load('t2.txt');  
%data_3 = load('L2_error_0.2.txt');  
%data_4 = load('L2_error_1.0.txt');  

% Extract the coordinate
x_coor = data_1(:,3);

% Extract the three types of solution
u_N_t1 = data_1(:,4);
u_t1   = data_1(:,5);
u_M_t1 = data_1(:,6);

u_N_t2 = data_2(:,4);
u_t2   = data_2(:,5);
u_M_t2 = data_2(:,6);


% Create figure with custom size
figure('Position', [100, 100, 1800, 800]);

% The first subplot
% Custom position: [left, bottom, width, height]
ax1 = axes('Position', [0.07, 0.15, 0.4, 0.78]); 
hold on;
plot(x_coor, u_N_t1, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Numerical');
plot(x_coor, u_t1,   'k--', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Analytical');
plot(x_coor, u_M_t1, 'b-.', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Modified');
set(gca, 'FontSize', 13);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 28);
ylabel('$u(x,t)$', 'Interpreter', 'latex', 'FontSize', 28);
title('Lax-Wendroff scheme solution at $t = 0.5$', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
grid on;

% The second subplot
ax2 = axes('Position', [0.52, 0.15, 0.4, 0.78]); 
hold on;
plot(x_coor, u_N_t2, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Numerical');
plot(x_coor, u_t2,   'k--', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Analytical');
plot(x_coor, u_M_t2, 'b-.', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Modified');
set(gca, 'FontSize', 13);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 28);
ylabel('$u(x,t)$', 'Interpreter', 'latex', 'FontSize', 28);
title('Lax-Wendroff scheme solution at $t = 1.5$', 'Interpreter', 'latex', 'FontSize', 22);
legend('Location', 'northeast', 'FontSize', 10, 'Interpreter', 'latex');
grid on;
hold off;

% save as pdf
saveas(gcf, 'solution_profile_N160.pdf');

