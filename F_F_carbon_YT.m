% 清理工作区和命令窗口
clear;
clc;

% --- 1. 读取 Excel 文件 ---
% 假设 carbon.xlsx 在 MATLAB 的当前工作目录下
filename = 'carbon_YD.xlsx';
% 读取所有数据，不带表头
data = readmatrix(filename);

% --- 2. 提取 YD 数据 ---
YD_data = data(1:249, :);

% --- 3. 按照 Depth 排序 YD 数据 ---
% Depth 是第二列 (data(:,2))
% 按照 Depth 从大到小（从浅到深）排序
% 注意：这里的排序逻辑是为了确保数据点是按照从浅到深排列的，
% 但Y轴的显示方向将由后续的plot设置决定。
[~, sort_idx] = sort(YD_data(:, 2), 'descend'); % 获取排序索引
sorted_YD_data = YD_data(sort_idx, :);

% 提取排序后的数据列
depth = sorted_YD_data(:, 2);
ccarb = sorted_YD_data(:, 3);
corg = sorted_YD_data(:, 4);
delta = sorted_YD_data(:, 5);

% --- 4. 绘制子图 ---
figure; % 创建一个新的图窗
% 设置图窗大小为 600x600 像素
% [left, bottom, width, height]
% left 和 bottom 定义图窗左下角在屏幕上的位置。
set(gcf, 'Position', [100, 100, 600, 600]); 

% 定义颜色
color1 = [53/255, 61/255, 38/255];      % 1元人民币绿色系
color2 = [114/255, 127/255, 101/255];   % 1元人民币绿色系
color3 = [189/255, 203/255, 177/255];   % 1元人民币绿色系

% 定义标记大小
marker_size = 4; % 调整标记大小

% 子图 1: Depth vs Ccarb
subplot(1, 3, 1); % 1行3列的第一个子图
plot(ccarb, depth, 'o-', 'Color', color1, ...
     'MarkerFaceColor', color1, 'MarkerSize', marker_size); % 实心圆并缩小
% 移除 set(gca, 'YDir', 'reverse'); 以恢复默认Y轴方向
xlabel('C_{carb}');
ylabel('Depth');
title('YD: Depth vs C_{carb}');
grid on;

% 子图 2: Depth vs Corg
subplot(1, 3, 2); % 1行3列的第二个子图
plot(corg, depth, 's-', 'Color', color2, ...
     'MarkerFaceColor', color2, 'MarkerSize', marker_size); % 实心方块并缩小
% 移除 set(gca, 'YDir', 'reverse');
xlabel('C_{org}');
ylabel('Depth');
title('YD: Depth vs C_{org}');
grid on;

% 子图 3: Depth vs Delta
subplot(1, 3, 3); % 1行3列的第三个子图
plot(delta, depth, '^-', 'Color', color3, ...
     'MarkerFaceColor', color3, 'MarkerSize', marker_size); % 实心三角形并缩小
% 移除 set(gca, 'YDir', 'reverse');
xlabel('Delta');
ylabel('Depth');
title('YD: Depth vs Delta');
grid on;

% 调整图窗布局，使其更美观
sgtitle('YD Carbon Data Analysis'); % 设置整个图窗的标题