      clc; clear; close all;

% load data
data_1 = load('solutionN200.txt');  
%data_2 = load('t2.txt');  
%data_3 = load('L2_error_0.2.txt');  
%data_4 = load('L2_error_1.0.txt');  

% Extract the coordinate
x_coor = data_1(:,1);

% Extract the three types of solution
u    = data_1(:,2);
u_M  = data_1(:,3);
u_N  = data_1(:,4);
u_AS = data_1(:,5);

%u_N_t2 = data_2(:,4);
%u_t2   = data_2(:,5);
%u_M_t2 = data_2(:,6);


% Create figure with custom size
figure('Position', [100, 100, 1200, 500]);

% The first subplot
% Custom position: [left, bottom, width, height]
ax1 = axes('Position', [0.08, 0.15, 0.4, 0.75]); 
hold on;
plot(x_coor, u    , 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'u');
plot(x_coor, u_M  , 'k--', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'u_M');
plot(x_coor, u_N  , 'b-.', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'u_N');
%plot(x_coor, u_AS , 'm*' , 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'u_AS');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$u(x,t)$', 'Interpreter', 'latex', 'FontSize', 24);
title('Solution at $t = 0.5$', 'Interpreter', 'latex', 'FontSize', 18);
legend('Location', 'best', 'FontSize', 10, 'Interpreter', 'latex');
grid on;

% The second subplot
%ax2 = axes('Position', [0.55, 0.15, 0.4, 0.75]); 
%hold on;
%plot(x_coor, u_N_t2, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Numerical');
%plot(x_coor, u_t2,   'k--', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Analytical');
%plot(x_coor, u_M_t2, 'b-.', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'Modified');
%xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24);
%ylabel('$u(x,t)$', 'Interpreter', 'latex', 'FontSize', 24);
%title('Solution at $t = 1.5$', 'Interpreter', 'latex', 'FontSize', 14);
%legend('Location', 'best', 'FontSize', 10, 'Interpreter', 'latex');
%grid on;
%hold off;

% save as pdf
saveas(gcf, 'solution_profile_N20.pdf');


