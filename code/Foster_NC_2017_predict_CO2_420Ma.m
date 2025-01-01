% 读取数据
filePath = "Foster_NC-2017_predict_CO2_420Ma.xlsx";
opts = detectImportOptions(filePath, 'Range', 'A3:F842'); % 从第三行到第842行读取数据
opts.VariableNames = {'Age_Ma', 'pCO2_probability_maximum', 'lw95', 'lw68', 'up68', 'up95'}; % 设置正确的列名
data = readtable(filePath, opts);

% 提取相关列
age = data.Age_Ma;  % 年代 (Ma)
pCO2_max = data.pCO2_probability_maximum;  % pCO2概率最大值
lower_95 = data.lw95;  % 95%置信区间下限
upper_95 = data.up95;  % 95%置信区间上限
lower_68 = data.lw68;  % 68%置信区间下限
upper_68 = data.up68;  % 68%置信区间上限

% 颜色的RGB值
color_line = [46/255, 159/255, 121/255];  % #2E9F79 (绿色)

% 绘制图表
figure;
hold on;
% 绘制95%置信区间 (浅灰色)
fill([age; flipud(age)], [lower_95; flipud(upper_95)], [0.8, 0.8, 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 绘制68%置信区间 (深灰色)
fill([age; flipud(age)], [lower_68; flipud(upper_68)], [0.5, 0.5, 0.5], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
% 绘制pCO2最大值曲线 (黑色)
plot(age, pCO2_max, 'k', 'LineWidth', 1.8);

% 添加横坐标线：644 和 1372
ymin = 644;
ymax = 1372;
line([age(1), age(end)], [ymin, ymin], 'Color', color_line, 'LineStyle', '--', 'LineWidth', 1.5); % 绘制y = 644的线
line([age(1), age(end)], [ymax, ymax], 'Color', color_line, 'LineStyle', '--', 'LineWidth', 1.5); % 绘制y = 1372的线


% 设置图形属性
xlabel('Age (Ma)');
ylabel('CO2 (ppm)');
title('CO2 Concentration Change Over Time');
grid on;

% 翻转横坐标轴
set(gca, 'XDir', 'reverse'); % 逆序显示横坐标

% 添加图例
legend({'95% CI', '68% CI', 'pCO2', 'Suitable range'}, ...
    'Location', 'Best', 'Box', 'off', 'EdgeColor', 'none');

hold off;
