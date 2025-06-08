% -----------------------------------------------------------------------
% Author: Your Name
% Date  : 2025-04-14 17:36:38
% File  : plot_scheme.m
% -----------------------------------------------------------------------
clc; clear; close all;

% 加载数据
data_1 = readmatrix('t1.txt', 'NumHeaderLines', 1);
data_2 = readmatrix('t2.txt', 'NumHeaderLines', 1);
data_3 = readmatrix('t3.txt', 'NumHeaderLines', 1);
data_4 = readmatrix('t4.txt', 'NumHeaderLines', 1);

% 提取坐标和数值解
x_coor = data_1(:, 2);
u_N1 = data_1(:, 3);   % Lax-Wendroff
u_N2 = data_1(:, 4);   % van Leer
u_N3 = data_1(:, 5);   % SUPERBEE
u_analytical = data_1(:, 6); % 解析解

% 创建自定义大小图像
figure('Position', [100, 100, 1400, 600], 'Color', 'w'); % 白色背景
hold on;
grid on;
set(gca, 'FontSize', 16, 'LineWidth', 1.5);

% ==================================================================
% 关键修改1：调整绘图顺序（最后绘制Lax-Wendroff以置顶）
% ==================================================================
plot(x_coor, u_analytical, 'k-', 'LineWidth', 3, 'DisplayName', 'Analytical Solution'); % 黑色实线解析解
plot(x_coor, u_N2, 'm--s', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'van Leer Limiter');
plot(x_coor, u_N3, 'b-.^', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'SUPERBEE Limiter');
plot(x_coor, u_N1, 'r-o', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Lax-Wendroff'); % 红色加粗实线

% ==================================================================
% 关键修改2：优化坐标轴和图例
% ==================================================================
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$u(x, t=0.5)$', 'Interpreter', 'latex', 'FontSize', 24);
title('Numerical Solutions at $t=0.5$', 'Interpreter', 'latex', 'FontSize', 26);

% 调整图例位置和样式
legend('Location', 'northeast', 'FontSize', 18, 'Interpreter', 'latex',...
       'Box', 'off', 'EdgeColor', 'none');

% 设置坐标轴范围（根据数据动态调整）
xlim([min(x_coor), max(x_coor)]);
ylim([-0.1 1.2]); % 留出边距避免标记被截断

% ==================================================================
% 关键修改3：增强网格和刻度可见性
% ==================================================================
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
box on;

% 保存为矢量图（避免像素化）
set(gcf, 'PaperPositionMode', 'auto');
print('Advection_Solutions.pdf', '-dpdf', '-r600');


