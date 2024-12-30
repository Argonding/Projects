% 读取Excel文件
filename = 'pbdb_biodiversity_data_from_540Ma.xlsx';
data = readtable(filename);

% 提取所需的列
max_ma = data.max_ma;        % 最大年龄（开始时间）
min_ma = data.min_ma;        % 最小年龄（结束时间）
sampled_in_bin = data.sampled_in_bin; % 物种多样性数据

% 创建一个新的时间轴，用于绘图
geological_time = (max_ma + min_ma) / 2; % 以每个区间的中点作为地质时间

% 只保留420 Ma到0 Ma之间的数据
filter = geological_time >= 0 & geological_time <= 420;
geological_time_filtered = geological_time(filter);
sampled_in_bin_filtered = sampled_in_bin(filter);

% 绘制曲线
figure;
plot(geological_time_filtered, sampled_in_bin_filtered, 'LineWidth', 2, 'Color', [0.545, 0, 0.071]); % 设置北大红颜色
xlabel('Geological Time (Ma)');
ylabel('Number of Genera');
title('Biodiversity Changes from 420 Ma to 0 Ma');

% 设定横轴范围从0到420，并设置逆序显示
set(gca, 'XDir', 'reverse'); % 逆序显示X轴
xlim([0 450]); % 横轴范围，从0到420
grid on;

% 设置字体
set(gca, 'FontSize', 12); % 设置字体大小
