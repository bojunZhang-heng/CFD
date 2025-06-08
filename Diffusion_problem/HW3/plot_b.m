clc; clear; close all;

% 加载数据
data_1 = load('L1_error_0.2.txt');  
data_2 = load('L1_error_1.0.txt');  
data_3 = load('L2_error_0.2.txt');  
data_4 = load('L2_error_1.0.txt');  

% 提取网格分辨率和误差
N = data_2(:,1);
L1_error_1 = data_1(:,3);
L1_error_2 = data_2(:,3);
L2_error_1 = data_3(:,3);
L2_error_2 = data_4(:,3);

% 创建图形窗口
figure('Position', [100, 100, 1000, 800]);

% 第一个子图
subplot(1,2,1);
cla reset; % 清除坐标轴
hold on;
h1 = loglog(N, L1_error_1, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'L1');
h2 = loglog(N, L2_error_1, 'b-s', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'L2');
xlabel('Grid Resolution (N)'); ylabel('Error');
title('$$\tau = 0.2$$', 'Interpreter', 'latex', 'FontSize', 12);
legend([h1, h2], 'Location', 'best');
grid on;
xticks([8 16 32 64 128]); 
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

% 第二个子图
subplot(1,2,2);
cla reset; % 清除坐标轴
hold on;
h3 = loglog(N, L1_error_2, 'r-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'L1');
h4 = loglog(N, L2_error_2, 'b-s', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'L2');
xlabel('Grid Resolution (N)'); ylabel('Error');
title('$$\tau = 1.0$$', 'Interpreter', 'latex', 'FontSize', 12);
legend([h3, h4], 'Location', 'best');
grid on;
xticks([8 16 32 64 128]); 
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

% 保存图形
saveas(gcf, 'error_profile.pdf');
