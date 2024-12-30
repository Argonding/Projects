% Step 1: 读取Excel数据
data = readtable('dataset_CO2_biodiv_from_420Ma.xlsx');

% Step 2: 根据地质时间从大到小排序
sortedData = sortrows(data, 'GeologicalTime_Ma', 'descend');

% Step 3: 绘制图表
figure;

% 绘制CO2预测及置信区间
yyaxis right;
% 绘制置信区间
hold on;
fill([sortedData.GeologicalTime_Ma; flip(sortedData.GeologicalTime_Ma)], ...
     [sortedData.lw95; flip(sortedData.up95)], [88/255, 182/255, 233/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
fill([sortedData.GeologicalTime_Ma; flip(sortedData.GeologicalTime_Ma)], ...
     [sortedData.lw68; flip(sortedData.up68)], [62/255, 144/255, 191/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
% 绘制CO2最大预测值线
plot(sortedData.GeologicalTime_Ma, sortedData.CO2_predict, '-', 'Color', [62/255, 79/255, 148/255], 'LineWidth', 1.5);

ylabel('CO2 Prediction(ppm)');
ax = gca;
ax.YColor = [62/255, 79/255, 148/255];  % 设置CO2纵轴为蓝色

% 然后绘制生物多样性曲线
yyaxis left;
plot(sortedData.GeologicalTime_Ma, sortedData.Smoothed_Biodiversity, '-', 'Color', [237/255, 141/255, 90/255], 'LineWidth', 1.5);
set(gca, 'XDir', 'reverse');  % 使x轴逆序
ylabel('Biodiversity(number of genera)');
ax.YColor = [237/255, 141/255, 90/255];  % 设置生物多样性纵轴为橙色

% 完成图表
xlabel('Geological Time (Ma)');
title('Biodiversity and CO2 Concentration over Geological Time', 'FontSize', 16); % 调整标题字体大小为16
grid on;

% 设置图例的位置，放置在上方正中间
legend('Biodiversity(number of genera)', '95% Confidence Interval', '68% Confidence Interval', 'CO2 Prediction(ppm)', ...
       'Location', 'north', 'Orientation', 'horizontal');

% 绘制适宜区间CO2浓度为644和1372的横线（加到右侧纵轴上）
yyaxis right;
h1 = yline(644, 'k--', 'LineWidth', 1.5);
h2 = yline(1372, 'k--', 'LineWidth', 1.5);

% 修改图例，只显示"Suitable CO2 Range"为标签
legend('Biodiversity(number of genera)', '95% Confidence Interval', '68% Confidence Interval', ...
       'CO2 Prediction(ppm)', 'Suitable CO2 Range', 'Location', 'north', 'Orientation', 'horizontal');