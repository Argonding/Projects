% 定义Excel文件的路径
filename = 'initial_LiXiang_dataset.xlsx';

% 读取Excel文件中的数据
data = readtable(filename);

% 提取'Surface_temperature'和'Precipitation'两列的数据
surfaceTemp = data.Surface_temperature;
CO2 = xlsread(filename, 'Sheet1', 'C2:C56'); % 读取第三列第2到56行的数据;

% 对CO2列取自然对数
log_CO2 = log(CO2);

% 将数据转换为表格形式，只包含用于回归的变量
tbl = table(surfaceTemp, 'VariableNames', {'Surface_temperature'});
tbl.log_CO2 = log_CO2;

% 进行线性回归分析
linearModel = fitlm(tbl, 'log_CO2 ~ Surface_temperature');

% 显示回归模型的统计信息
disp(linearModel)

% 获取线性回归的系数
coefficients = linearModel.Coefficients.Estimate;
slope = coefficients(2); % 斜率
intercept = coefficients(1); % 截距

% 绘制散点图和回归线
figure;
% 设置y轴为对数坐标
scatter(surfaceTemp, CO2, 'filled');
hold on;
set(gca,'YScale','log'); % 设置y轴为对数坐标

% 计算回归线上的预测值
surfaceTempFit = linspace(min(surfaceTemp), max(surfaceTemp), 100).';
tblFit = table(surfaceTempFit, 'VariableNames', {'Surface_temperature'});

% 使用表格形式的数据来预测
log_CO2Fit = predict(linearModel, tblFit);

% 回归线应基于log_CO2Fit的预测值
CO2Fit = exp(log_CO2Fit); % 反变换回原始尺度

% 绘制回归线
plot(surfaceTempFit, CO2Fit, '-r');

% 添加标题和标签
xlabel('Surface Temperature(℃)');
ylabel('CO2(/280ppmv)');
title('Linear Regression of Surface Temperature vs. log(CO2(/280ppmv))');

% 在图中标注线性回归方程
annotation('textbox', [0.2 0.7 0.1 0.1], ...
    'String', sprintf('log(y) = %.2fx + %.2f', slope, intercept), ...
    'Interpreter', 'latex', 'FontSize', 15);

legend('Data', 'Linear Fit in Log Scale');
hold off;