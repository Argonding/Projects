% 参数设置
Um = 2.0; % 浮游植物最大营养盐摄入率
kN = 1.0; % 浮游植物吸收DIN的半饱和常数
Gm = 0.2; % 浮游动物最大摄食率
lambda = 0.2; % Ivlev摄入常数
gamma = 0.2; % 浮游动物生长系数
theta = 0.4; % 浮游动物排泄系数
Is = 82.0; % 平均海水表面光辐射强度
I0 = 150.0; % 最优光强
Qg10 = 2.08; % 生长温度依赖系数
Qh10 = 3.1; % 摄食温度依赖系数
Mp = 0.05; % 浮游植物基础死亡率
MZ = 0.25; % 浮游动物基础死亡率
e = 0.05; % 碎屑再矿化率
k_H = 3.3e-2; % CO2的亨利常数，单位mol/L·atm
P_total = 1; % 总大气压，单位atm
k_CO2 = 0.00005; % CO2的半饱和常数
CO2_max = 0.001; % CO2浓度的抑制阈值
alpha = 0.000005; % 酸化程度影响系数

% 读取CO2大气浓度数据
filename = "Foster_NC-2017_predict_CO2_420Ma.xlsx";
CO2_data = xlsread(filename, 'LOESS Fit', 'B3:B842'); % 读取第二列第3到842行的数据:pCO2_probability_maximum

% 初始条件
N0 = 0.99; % 营养盐初始浓度
P0 = 0.01; % 浮游植物初始浓度
Z0 = 0.001; % 浮游动物初始浓度
D0 = 0.99; % 碎屑初始浓度
tspan = [0 200]; % 模拟时长从0到200天
y0 = [N0, P0, Z0, D0]; % 初始状态向量

% 存储崩溃和未崩溃的CO2浓度
crash_CO2 = [];
no_crash_CO2 = [];

for i = 1:30
    CO2_air_ppm = CO2_data(i); 
    
    % 计算大气中CO2的分压
    P_CO2 = CO2_air_ppm * P_total / 1e6;
    
    % 计算CO2在水中的浓度
    CO2_sea = k_H * P_CO2;
    
    % 计算酸化程度
    acidification = P_CO2;
    
    % 定义微分方程
    odefun_with_CO2 = @(t, y) npzd_ode_with_CO2(t, y, Um, kN, Gm, lambda, gamma, theta, Is, I0, Qg10, Qh10, Mp, MZ, e, CO2_sea, k_CO2, CO2_max, CO2_air_ppm, alpha, acidification);
    
    % 使用ode15s求解
    [t, y] = ode15s(odefun_with_CO2, tspan, y0);
    
    % 取模拟末尾部分的浮游植物浓度
    P_end = y(end, 2); % 最后一个时间点的浮游植物浓度
    P_initial = P0;    % 初始的浮游植物浓度
    
    % 获取最后50个时间点的浮游植物浓度
    if length(y) >= 50
        P_last_section = y(end-49:end, 2);
    else
        P_last_section = y(:, 2);
    end

    % 判断浮游植物是否崩溃：如果最后50个时间点的最大值小于初始值的30%
    if max(P_last_section) <= 0.3 * P_initial
        crash_CO2 = [crash_CO2, CO2_data(i)]; % 记录崩溃时的CO2浓度
    else
        no_crash_CO2 = [no_crash_CO2, CO2_data(i)]; % 记录未崩溃的CO2浓度
    end
    
end

for i = 31:length(CO2_data)
    CO2_air_ppm = CO2_data(i); 
    
    % 计算大气中CO2的分压
    P_CO2 = CO2_air_ppm * P_total / 1e6;
    
    % 计算CO2在水中的浓度
    CO2_sea = k_H * P_CO2;
    
    % 计算酸化程度
    acidification = P_CO2;
    
    % 定义微分方程
    odefun_with_CO2 = @(t, y) npzd_ode_with_CO2(t, y, Um, kN, Gm, lambda, gamma, theta, Is, I0, Qg10, Qh10, Mp, MZ, e, CO2_sea, k_CO2, CO2_max, CO2_air_ppm, alpha, acidification);
    
    % 使用ode15s求解
    [t, y] = ode15s(odefun_with_CO2, tspan, y0);
    
    % 取模拟末尾部分的浮游植物浓度
    P_end = y(end, 2); % 最后一个时间点的浮游植物浓度
    P_initial = P0;    % 初始的浮游植物浓度
    
    % 获取最后50个时间点的浮游植物浓度
    if length(y) >= 50
        P_last_section = y(end-49:end, 2);
    else
        P_last_section = y(:, 2);
    end

    % 判断浮游植物是否稳定并检查是否达到崩溃条件
    if P_end <= 0.3 * P_initial
        crash_CO2 = [crash_CO2, CO2_data(i)]; % 记录崩溃时的CO2浓度
    else
        no_crash_CO2 = [no_crash_CO2, CO2_data(i)]; % 记录未崩溃的CO2浓度
    end
    
end

% 排序未崩溃的CO2浓度（按从高到低排序）
sorted_no_crash_CO2 = sort(no_crash_CO2, 'descend'); % 按从高到低排序

% 排序崩溃的CO2浓度（按从高到低排序）
sorted_crash_CO2 = sort(crash_CO2, 'descend'); % 按从高到低排序

% 计算第一个和第三个四分位数（Q1和Q3）
Q1 = prctile(sorted_no_crash_CO2, 25);
Q3 = prctile(sorted_no_crash_CO2, 75);

% 计算四分位距（IQR）
IQR = Q3 - Q1;

% 定义非离群点数据的下边界和上边界
lower_bound = Q1 - 1.5 * IQR;
upper_bound = Q3 + 1.5 * IQR;

% 过滤数据以排除离群点
filtered_no_crash_CO2 = sorted_no_crash_CO2(sorted_no_crash_CO2 >= lower_bound & sorted_no_crash_CO2 <= upper_bound);

% Step 1: 读取Excel数据
data = readtable('dataset_CO2_biodiv_from_420Ma.xlsx');

% Step 2: 根据地质时间从大到小排序
sortedData = sortrows(data, 'GeologicalTime_Ma', 'descend');

% Step 3: 绘制原始图表
figure;

% 绘制生物多样性曲线 (移到右侧)
yyaxis right;
biodiversityCurve = plot(sortedData.GeologicalTime_Ma, sortedData.Smoothed_Biodiversity, '-', ...
     'Color', [237/255, 141/255, 90/255], 'LineWidth', 1.5);
set(gca, 'XDir', 'reverse');  % 使x轴逆序
ylabel('Biodiversity (number of genera)');
ax = gca;
ax.YColor = [237/255, 141/255, 90/255];  % 设置生物多样性纵轴为橙色

% 绘制CO2预测及置信区间 (移到左侧)
yyaxis left;
hold on;
fill95 = fill([sortedData.GeologicalTime_Ma; flip(sortedData.GeologicalTime_Ma)], ...
     [sortedData.lw95; flip(sortedData.up95)], [88/255, 182/255, 233/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
fill68 = fill([sortedData.GeologicalTime_Ma; flip(sortedData.GeologicalTime_Ma)], ...
     [sortedData.lw68; flip(sortedData.up68)], [62/255, 144/255, 191/255], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
co2PredictionLine = plot(sortedData.GeologicalTime_Ma, sortedData.CO2_predict, '-', ...
     'Color', [62/255, 79/255, 148/255], 'LineWidth', 1.5);
ylabel('CO2 Prediction (ppm)');
ax = gca;
ax.YColor = [62/255, 79/255, 148/255];  % 设置CO2纵轴为蓝色

% 绘制适宜区间CO2浓度为644和1372的横线（加到左侧纵轴上）
line1 = yline(644, 'Color', [0.07, 0.35, 0.40], 'LineStyle', '--', 'LineWidth', 2.4); % 改为#845ec2
line2 = yline(1372, 'Color', [0.07, 0.35, 0.40], 'LineStyle', '--', 'LineWidth', 2.4); % 改为#845ec2

% 添加基于筛选后CO2数据最大值和最小值的新横线
if ~isempty(filtered_no_crash_CO2)
    line3 = yline(max(filtered_no_crash_CO2), 'Color', [0.60, 0.24, 0.18], 'LineStyle', '--', 'LineWidth', 2.4); % 新的最大值线，颜色#008c81
    line4 = yline(min(filtered_no_crash_CO2), 'Color', [0.60, 0.24, 0.18], 'LineStyle', '--', 'LineWidth', 2.4); % 新的最小值线，颜色#008c81
end

% Step 4 & 5: 在x = 425Ma 和 x = 435Ma 处绘制一个垂直于x轴的直线，并投点（统一崩溃和未崩溃CO2数据的颜色）

% 绘制x=425Ma处的垂直线及投点
x1 = 425;
vline1 = xline(x1, '--');
hold on;

% 绘制未崩溃的CO2数据，使用新颜色表示
if ~isempty(filtered_no_crash_CO2)
    nonCrashPoints = plot(repmat(x1, length(filtered_no_crash_CO2), 1), filtered_no_crash_CO2, 'o', ...
         'MarkerFaceColor', [62/255, 79/255, 148/255], 'MarkerEdgeColor', [62/255, 79/255, 148/255]); % 新颜色
end

% 绘制x=435Ma处的垂直线及投点
x2 = 435;
vline2 = xline(x2, '--');

% 分别绘制崩溃和未崩溃的数据点，保持与x=425Ma一致的颜色方案
if ~isempty(sorted_crash_CO2)
    crashPoints = plot(repmat(x2, length(sorted_crash_CO2), 1), sorted_crash_CO2, 'o', ...
         'MarkerFaceColor', [237/255, 141/255, 90/255], 'MarkerEdgeColor', [237/255, 141/255, 90/255]); % 新颜色
end
if ~isempty(sorted_no_crash_CO2)
    plot(repmat(x2, length(sorted_no_crash_CO2), 1), sorted_no_crash_CO2, 'o', ...
         'MarkerFaceColor', [62/255, 79/255, 148/255], 'MarkerEdgeColor', [62/255, 79/255, 148/255]); % 新颜色
end

% 确保所有绘图命令后关闭hold
hold off;

% 添加标题
title('Biodiversity v.s. CO2 in the atmosphere from 420Ma', 'FontSize', 14);

% 添加图例
legend([biodiversityCurve, co2PredictionLine, fill95, fill68, nonCrashPoints, crashPoints, line1, line3], ...
       {'Biodiversity', 'CO2 Prediction(probability maximum)', '95%CI CO2', '68%CI CO2', 'Non-Crash CO2 Data', 'Crash CO2 Data', 'Suitable CO2 Range 1: 644 - 1372 ppm(Database A)', ...
        ['Suitable CO2 Range 2: ' num2str(min(filtered_no_crash_CO2)) ' - ' num2str(max(filtered_no_crash_CO2)) ' ppm(Database B)'], ...
       }, ...
       'Location', 'best');




% 定义微分方程
function dydt = npzd_ode_with_CO2(t, y, Um, kN, Gm, lambda, gamma, theta, Is, I0, Qg10, Qh10, Mp, MZ, e, CO2_concentration, k_CO2, C_max, CO2_air_ppm, alpha, acidification)
    N = y(1);
    P = y(2);
    Z = y(3);
    D = y(4);

    % 计算温度 T
    T = (log(CO2_air_ppm / 280) + 1.83) / 0.19;
    
    gT = Qg10^(T - 10); % 生长函数
    hT = Qh10^(T - 10); % 摄食函数
    
    % 光限制函数
    S = 10; 
    H = 3 * S;
    I_s = H * 1.51 * (1 - exp(-1.51 * H / S));
    fI = 1 / (1 - exp(-4.53)) * (1 - (1 / 4.53) * (I_s / I0) * (1 - exp(-4.53))) * (I_s / I0);
   
    % 调整CO2影响浮游植物生长的函数
    J_CO2 = (CO2_concentration / (k_CO2 + CO2_concentration)) * (1 - (CO2_concentration / C_max));
    
    % 考虑酸化程度对浮游植物和浮游动物死亡率的影响
    mortality_p = Mp * (1 + alpha * acidification);
    mortality_z = MZ * (1 + alpha * acidification);
    
    % 计算变化率
    dNdt = -Um * fI * gT * J_CO2 * P * (N / (kN + N)) + theta * Gm * hT * Z * (1 - exp(-lambda * P)) + e * D;
    dPdt = Um * fI * gT * J_CO2 * P * (N / (kN + N)) - Gm * hT * Z * (1 - exp(-lambda * P)) - mortality_p * P;
    dZdt = gamma * Gm * hT * Z * (1 - exp(-lambda * P)) - mortality_z * Z;
    dDdt = (1 - gamma - theta) * Gm * hT * Z * (1 - exp(-lambda * P)) + mortality_p * P + mortality_z * Z - e * D;
    
    dydt = [dNdt; dPdt; dZdt; dDdt];
end


