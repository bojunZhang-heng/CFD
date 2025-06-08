clc; clear; close all;

%------------------------------------------------------------------------
% Load data
data_1 = load('velocity_1.txt');  
data_2 = load('velocity_2.txt');  
data_3 = load('velocity_3.txt');  
data_4 = load('velocity_4.txt');  

% Extract dimensionless time 
time_1 = data_1(:,2);
time_2 = data_2(:,2);
time_3 = data_3(:,2);
time_4 = data_4(:,2);

% Extract velocity
velocity_1 = data_1(:,3);
velocity_2 = data_2(:,3);
velocity_3 = data_3(:,3);
velocity_4 = data_4(:,3);

%------------------------------------------------------------------------
% Plot in separate subplots
figure('Position', [100, 100, 1000, 800]); % Adjust figure size

% First plot
subplot(2,2,1);
plot(time_1, velocity_1, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 1);
xlabel('Time'); ylabel('Velocity');
title('CFL = 0.500');
grid on;

% Second plot
subplot(2,2,2);
plot(time_2, velocity_2, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 1);
xlabel('Time'); ylabel('Velocity');
title('CFL = 0.505');
grid on;

% Third plot
subplot(2,2,3);
plot(time_3, velocity_3, 'g--^', 'LineWidth', 1.5, 'MarkerSize', 1);
xlabel('Time'); ylabel('Velocity');
title('CFL = 0.510');
grid on;

% Fourth plot
subplot(2,2,4);
plot(time_4, velocity_4, 'm--d', 'LineWidth', 1.5, 'MarkerSize', 1);
xlabel('Time'); ylabel('Velocity');
title('CFL = 0.520');
grid on;

% Save as PDF
saveas(gcf, 'four_velocity_profiles.pdf');  







