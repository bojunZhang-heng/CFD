    clc; clear; close all;

% Open the file
filename = 'shear_stress_result_N32.txt';                   % Ensure the file is in the working directory
fid = fopen(filename, 'r');

if fid == -1
    error('File not found or could not be opened.');
end

% Read the entire file into a cell array
data = textscan(fid, '%f');
fclose(fid);

% Extract numeric values
tau = data{1};    % Forth column 

% Determine the number of data points per block
N = numel(tau) / 2;

if mod(numel(tau), 2) ~= 0
    error('The number of data lines is not even. Check the file format.');
end

% Split the data into two sets
tt = [0.2; 1.0];

tau_analytical = tau(1:N);            % First set of analytical solutions
tau_numerical  = tau(N+1:end);        % Second set of analytical solutions


%------------------------------------------------------------------------
% Plot
figure('Position', [100, 100, 900, 700]); % Larger figure size
hold on;

% Numerical solutions (markers)
plot(tt, tau_analytical, '-o', 'Color', [0 0.447 0.741], 'MarkerFaceColor', [0 0.447 0.741], ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', '$Analytical$');
plot(tt, tau_numerical,  '-s', 'Color', [0.85 0.325 0.098], 'MarkerFaceColor', [0.85 0.325 0.098], ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', '$Numerical$');

% Labels and Title
xlabel('dimenionless time', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('tau', 'FontSize', 14, 'Interpreter', 'latex');
title('shear stress evolution', 'FontSize', 16, 'Interpreter', 'latex');

% Legend
legend('show', 'Location', 'northeastoutside', 'FontSize', 12, 'Interpreter', 'latex');

% Grid
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

hold off;

% Save as PDF
saveas(gcf, 'shear_stress_time_evolution_N32.pdf');  

