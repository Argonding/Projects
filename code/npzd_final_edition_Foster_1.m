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

% 创建新的图形窗口
figure;

% 绘制过滤后的未崩溃CO2数据，使用新的颜色表示
hold on; % 保持当前图形以便绘制多个数据集
plot(filtered_no_crash_CO2, zeros(size(filtered_no_crash_CO2)), 'o', 'MarkerFaceColor', [62/255, 79/255, 148/255], 'MarkerEdgeColor', [62/255, 79/255, 148/255]); % 蓝色圆圈

grid on

% 设置X轴标签
xlabel('CO2浓度（ppm）');
% 移除Y轴刻度，因为我们只关心X轴上的分布
set(gca, 'YTick', []);
% 添加标题
title('一维轴上的过滤后的未崩溃CO2浓度');
% 确保所有点都能显示在图中
axis tight;


% 创建新的图形窗口
figure;

% 绘制未崩溃的CO2数据，使用新的颜色表示
hold on; % 保持当前图形以便绘制多个数据集
plot(sorted_no_crash_CO2, zeros(size(sorted_no_crash_CO2)), 'o', 'MarkerFaceColor', [62/255, 79/255, 148/255], 'MarkerEdgeColor', [62/255, 79/255, 148/255]); % 蓝色圆圈
% 绘制崩溃的CO2数据，使用新的颜色表示
plot(sorted_crash_CO2, zeros(size(sorted_crash_CO2)), 'o', 'MarkerFaceColor', [237/255, 141/255, 90/255], 'MarkerEdgeColor', [237/255, 141/255, 90/255]); % 红色圆圈

grid on

% 设置X轴标签
xlabel('CO2 Concentration (ppm)');
% 移除Y轴刻度，因为我们只关心X轴上的分布
set(gca, 'YTick', []);
% 添加标题
title('一维轴上的崩溃/未崩溃CO2浓度');
% 添加图例
legend('未崩溃CO2浓度', '崩溃CO2浓度');
% 确保所有点都能显示在图中
axis tight;

% 创建新的图形窗口
figure; % 创建一个新的图形窗口

% 为未崩溃的CO2浓度绘制箱线图和数据点
subplot(1,2,1); % 在一个1x2的网格中创建第一个子图
boxplot(sorted_no_crash_CO2, 'Colors', [62/255, 79/255, 148/255], 'Widths', 0.5); % 去掉LineWidth参数
hold on; % 保持当前图形
scatter(repmat(1, size(sorted_no_crash_CO2)), sorted_no_crash_CO2, 'filled', 'MarkerFaceColor', [62/255, 79/255, 148/255], 'MarkerEdgeColor', [62/255, 79/255, 148/255]); % 绘制数据点
title('未崩溃情况下的CO2浓度分布');
xlabel(' ');
ylabel('大气CO2浓度');
ylim([0 2200]); % 设置纵坐标范围

% 为崩溃的CO2浓度绘制箱线图和数据点
subplot(1,2,2); % 在同一个1x2的网格中创建第二个子图
boxplot(sorted_crash_CO2, 'Colors', [237/255, 141/255, 90/255], 'Widths', 0.5); % 去掉LineWidth参数
hold on; % 保持当前图形
scatter(repmat(1, size(sorted_crash_CO2)), sorted_crash_CO2, 'filled', 'MarkerFaceColor', [237/255, 141/255, 90/255], 'MarkerEdgeColor', [237/255, 141/255, 90/255]); % 绘制数据点
title('崩溃情况下的CO2浓度分布');
xlabel(' ');
ylabel('大气CO2浓度');
ylim([0 2200]); % 设置纵坐标范围


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


