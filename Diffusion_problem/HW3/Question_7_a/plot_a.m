    clc; clear; close all;

% Load data
data = load('relative_error.txt');  

% Extract columns
y = data(:,1);  
time = data(:,2);
errors = data(:,5);
tt = [];
trunc_error_y0 = [];
trunc_error_y1 = [];
whole_error_y0 = [];
whole_error_y1 = [];
tt = [];

%------------------------------------------------------------------------
% Find indices for y = 0 and y = -0.5
idx_y0 = (y == 0);
idx_y_neg_half = (y == -0.5);

%------------------------------------------------------------------------
% Extract dimensionless time
count_0 = 0;                                 % Initialize counter
for ii = 1 : length(time)
    if mod(ii, 4) == 0
      count_0 = count_0 + 1;
      tt(count_0,1) = time(ii);
    end
end

%------------------------------------------------------------------------
% Extract error
count_1 = 0;
count_2 = 0;
count_3 = 0;
count_4 = 0;
for ii = 1 : length(errors)
  if mod(ii,4) == 1
    count_1 = count_1 + 1;
    trunc_error_y0(count_1,1) = errors(ii);
  end

  if mod(ii,4) == 2
    count_2 = count_2 + 1;
    whole_error_y0(count_2,1) = errors(ii);
  end

  if mod(ii,4) == 3
    count_3 = count_3 + 1;
    trunc_error_y1(count_3,1) = errors(ii);
  end


  if mod(ii,4) == 0
    count_4 = count_4 + 1;
    whole_error_y1(count_4,1) = errors(ii);
  end
end
%------------------------------------------------------------------------
% Plot
figure('Position', [100, 100, 900, 700]); % Larger figure size
hold on;

% Numerical solutions (markers)
plot(tt, trunc_error_y0, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 1, 'DisplayName', 'Trunc\_Error at y=0');
plot(tt, trunc_error_y1, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 1, 'DisplayName', 'Trunc\_Error at y=-0.5L');
plot(tt, whole_error_y0, 'g--^', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', 'Whole\_Error at y=0');
plot(tt, whole_error_y1, 'm--d', 'LineWidth', 1.5, 'MarkerSize', 1, 'DisplayName', 'Whole\_Error at y=-0.5L');

% Customize plot
xlabel('Dimensionless Time', 'FontSize', 12);
ylabel('Error', 'FontSize', 12);
title('Error Evolution in Planar Poiseuille Flow', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
grid on;
box on;
set(gca, 'FontSize', 10);

% Adjust axis limits
xlim([min(tt)-0.01, max(tt)+0.01]);
ylim([min([trunc_error_y0; trunc_error_y1; whole_error_y0; whole_error_y1])*0.9, ...
      max([trunc_error_y0; trunc_error_y1; whole_error_y0; whole_error_y1])*1.1]);

hold off;

% Save as PDF
saveas(gcf, 'relative_error.pdf');  




