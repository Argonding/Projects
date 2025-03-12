% 定义常数
K1 = 10^-6.35; % 碳酸的第一级解离常数
K2 = 10^-10.33; % 碳酸的第二级解离常数
Kw = 10^-14; % 水的离子积常数
KCO2 = 10^-1.46; % CO2的溶解度常数
pCO2 = 10^-3.5; % 大气中CO2的分压

% pH范围
pH_range = 0:0.1:14;

% 初始化数组存储结果
log_M_H2CO3 = zeros(size(pH_range));
log_M_HCO3_minus = zeros(size(pH_range));
log_M_CO3_double_minus = zeros(size(pH_range));
log_a_H_plus = zeros(size(pH_range));
log_a_OH_minus = zeros(size(pH_range));
log_CT = zeros(size(pH_range));

% 计算不同pH值下的各组分活度
for i = 1:length(pH_range)
    pH = pH_range(i);
    a_H_plus = 10^-pH; % 氢离子活度
    
    % 计算碳酸、碳酸氢根和碳酸根的活度
    M_H2CO3 = KCO2 * pCO2;
    log_M_H2CO3(i) = log10(M_H2CO3);
    
    M_HCO3 = (K1 * M_H2CO3) / a_H_plus;
    log_M_HCO3_minus(i) = log10(M_HCO3);
    
    M_CO3 = (K1 * K2 * M_H2CO3) / (a_H_plus^2);
    log_M_CO3_double_minus(i) = log10(M_CO3);
    
    % 氢离子和氢氧根活度
    log_a_H_plus(i) = -pH;
    log_a_OH_minus(i) = -(14 - pH); % 利用水的离子积计算OH-浓度
    
    % 总无机碳浓度
    CT = M_H2CO3 + M_HCO3 + M_CO3;
    log_CT(i) = log10(CT);
end

% 绘制图形
figure;
plot(pH_range, log_M_H2CO3, 'Color', [35/255 51/255 66/255], 'LineWidth', 2, 'DisplayName', 'H_2CO_3');
hold on;
plot(pH_range, log_M_HCO3_minus, 'Color', [74/255 90/255 105/255], 'LineWidth', 2, 'DisplayName', 'HCO_3^-');
plot(pH_range, log_M_CO3_double_minus, 'Color', [167/255 174/255 190/255], 'LineWidth', 2, 'DisplayName', 'CO_3^{2-}');
plot(pH_range, log_a_H_plus, 'Color', [184/255 164/255 157/255], 'LineWidth', 2, 'DisplayName', 'H^+');
plot(pH_range, log_a_OH_minus, 'Color', [181/255 179/255 166/255], 'LineWidth', 2, 'DisplayName', 'OH^-');
plot(pH_range, log_CT, 'Color', [183/255 140/255 124/255], 'LineWidth', 2, 'DisplayName', 'Total DIC');

xlabel('pH');
ylabel('Log of Activity');
legend('show');
title('Activity of Carbon Species vs pH', 'FontSize',20);
grid on;

% 设置y轴的显示范围为从-8到0
ylim([-8 0]);
