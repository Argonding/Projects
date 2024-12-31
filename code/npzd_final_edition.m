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
filename = "initial_LiXiang_dataset.xlsx";
CO2_data = xlsread(filename, 'Sheet1', 'C2:C56'); % 读取第三列第2到56行的数据

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

% 创建第一张图
figure;
set(gcf, 'Position', [100, 100, 1400, 1000]); % 设置窗口大小
for i = 1:30
    CO2_air_ppm = CO2_data(i) * 280; % 将CO2浓度×280代入模型
    
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
    P_last_section = y(end-49:end, 2);
    
    % 判断浮游植物是否崩溃：如果最后50个时间点的最大值小于初始值的30%
    if max(P_last_section) <= 0.3 * P_initial
        crash_CO2 = [crash_CO2, CO2_data(i)]; % 记录崩溃时的CO2浓度
    else
        no_crash_CO2 = [no_crash_CO2, CO2_data(i)]; % 记录未崩溃的CO2浓度
    end
    
    % 创建子图：6行5列的布局，增加宽度
    subplot(6, 5, i); % 6行5列的子图布局
    plot(t, y);
    xlabel('Time (days)');
    ylabel('Concentration');
    title(['CO2: ', num2str(CO2_data(i)), ' ppm']);
    
    % 设置Y轴范围
    ylim([0, 1]);
    

end
sgtitle('NPZD Model Simulation for CO2(×280) Concentrations');

% 创建第二张图
figure;
set(gcf, 'Position', [100, 100, 1400, 1000]); % 设置窗口大小
for i = 31:length(CO2_data)
    CO2_air_ppm = CO2_data(i) * 280; % 将CO2浓度×280代入模型
    
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
    
    % 判断浮游植物是否稳定并检查是否达到崩溃条件
    if P_end <= 0.3 * P_initial
        crash_CO2 = [crash_CO2, CO2_data(i)]; % 记录崩溃时的CO2浓度
    else
        no_crash_CO2 = [no_crash_CO2, CO2_data(i)]; % 记录未崩溃的CO2浓度
    end
    
    % 创建子图：5行5列的布局
    subplot(5, 5, i-30); % 5行5列的子图布局
    plot(t, y);
    xlabel('Time (days)');
    ylabel('Concentration');
    title(['CO2: ', num2str(CO2_data(i)), ' ppm']);
    
    % 设置Y轴范围
    ylim([0, 1]);
    
    % 添加标签和图例
    if i == 31
        legend('Nutrient (N)', 'Phytoplankton (P)', 'Zooplankton (Z)', 'Detritus (D)', 'Location', 'best');
    end
end
sgtitle('NPZD Model Simulation for CO2(×280) Concentrations');

% 排序并输出崩溃的CO2浓度（按从高到低排序）
sorted_crash_CO2 = sort(crash_CO2, 'descend'); % 按从高到低排序
disp('按从高到低顺序输出的导致系统崩溃的CO2浓度值：');
disp(sorted_crash_CO2);

% 排序并输出未崩溃的CO2浓度（按从高到低排序）
sorted_no_crash_CO2 = sort(no_crash_CO2, 'descend'); % 按从高到低排序
disp('按从高到低顺序输出的未导致系统崩溃的CO2浓度值：');
disp(sorted_no_crash_CO2);


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

