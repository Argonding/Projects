% 读取数据
file_path = 'C:\Users\Argon\Desktop\final_project\code_and_dataset\dataset_CO2_biodiv_from_420Ma.xlsx';
data = readtable(file_path);

% 提取地质时间、Smoothed_Biodiversity 和 CO2_predict 三列
time = data.GeologicalTime_Ma;
biodiversity = data.Smoothed_Biodiversity;
CO2 = data.CO2_predict;

% 筛选出中生代（251 Ma - 66 Ma）的数据
mesozoic_mask = (time >= 66) & (time <= 251);
mesozoic_time = time(mesozoic_mask);
mesozoic_biodiversity = biodiversity(mesozoic_mask);
mesozoic_CO2 = CO2(mesozoic_mask);

% 设置滞后期范围：从0.5 Ma到10 Ma，步长为0.5 Ma
lags = 0.5:0.5:10;  % 滞后期范围
lag_correlations = zeros(length(lags), 1); % 存储每个滞后期的相关性

% 计算不同滞后期的相关性
for i = 1:length(lags)
    lag = lags(i); 
    
    % 滞后生物多样性数据，直接向后平移
    shifted_biodiversity = mesozoic_biodiversity(round(lag*2)+1:end); % 滞后后的生物多样性数据（乘2是为了转换为索引）
    
    % CO2数据与滞后生物多样性数据对齐
    shifted_CO2 = mesozoic_CO2(1:end-round(lag*2)); % 取滞后后的CO2数据（乘2是为了转换为索引）
    
    % 计算滞后期的皮尔逊相关系数
    lag_correlations(i) = corr(shifted_biodiversity, shifted_CO2, 'Type', 'Pearson');
end

% 找到最小的相关系数，表示最强的负相关性
[~, min_index] = min(lag_correlations); % 找到最小的相关系数的索引
best_lag = lags(min_index);  % 对应的滞后期

% 输出最佳滞后期
disp(['最佳滞后期为 ', num2str(best_lag), ' Ma, 对应的皮尔逊相关系数为 ', num2str(lag_correlations(min_index))]);


% 绘制滞后期与皮尔逊相关系数的图像
figure;
plot(lags, lag_correlations, '-o');
xlabel('滞后期（Ma）');
ylabel('皮尔逊相关系数');
title('滞后期与皮尔逊相关系数的关系');
grid on;
