% 碳酸的第一和第二解离常数
K1 = 10^(-6.35); % 第一解离常数
K2 = 10^(-10.33); % 第二解离常数

% 总无机碳浓度 c (mol/L)
c = 2; 

% pH范围
pH_values = linspace(0, 14, 280); % 从0到14，共280个点

% 初始化数组保存结果
a_H2CO3 = zeros(size(pH_values));
a_HCO3 = zeros(size(pH_values));
a_CO3 = zeros(size(pH_values));

for i = 1:length(pH_values)
    pH = pH_values(i);
    H_concentration = 10^(-pH); % 氢离子浓度
    
    % 统一分母 D
    D = 10^(-2*pH) + K1 * H_concentration + K1*K2;
    
    % 各物种活度计算
    a_H2CO3(i) = c * 10^(-2*pH) / D;
    a_HCO3(i) = c * K1 * H_concentration / D;
    a_CO3(i) = c * K1 * K2 / D;
end

% 创建一个新的图形窗口
figure;

plot(pH_values, a_H2CO3, 'Color', [53/255 61/255 38/255], 'LineWidth', 2, 'DisplayName', 'H_2CO_3');
hold on;
plot(pH_values, a_HCO3, 'Color', [114/255 127/255 101/255], 'LineWidth', 2, 'DisplayName', 'HCO_3^-');
plot(pH_values, a_CO3, 'Color', [189/255 203/255 177/255], 'LineWidth', 2, 'DisplayName', 'CO_3^{2-}');

% 添加图例、标题和轴标签
xlabel('pH');
ylabel('activity(mM)');
title('the activity of different inorganic carbon species in solution with respect to pH');
legend show;

% 显示网格线
grid on;

% 结束保持当前图形的状态
hold off;
