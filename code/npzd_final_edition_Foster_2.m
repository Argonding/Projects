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

% 读取真实温度和CO2浓度数据
lixFile = 'initial_LiXiang_dataset.xlsx';
lixData = readtable(lixFile);

fosterFile = 'Foster_NC-2017_predict_CO2_420Ma.xlsx';
sheets = sheetnames(fosterFile);
fosterData = [];

for i = 1:length(sheets)
    tempData = readtable(fosterFile, 'Sheet', sheets{i});
    fosterData = [fosterData; tempData];
end

timePoints = (410:-10:0)'; % 创建时间点列向量
realTemperatures = zeros(length(timePoints), 1); % 真实温度向量
realCO2Concentrations = zeros(length(timePoints), 1); % 真实CO2浓度向量

% 提取真实温度和CO2浓度值
for i = 1:length(timePoints)
    age = timePoints(i);
    idx = find(abs(fosterData.Age_Ma - age) == min(abs(fosterData.Age_Ma - age)), 1);
    
    realCO2Concentrations(i) = fosterData.pCO2_probability_maximum(idx);
    realTemperatures(i) = lixData.Surface_temperature(lixData.Year_Ma == age);
end

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

% 模拟循环
for i = 1:length(realCO2Concentrations)
    CO2_air_ppm = realCO2Concentrations(i);
    T = realTemperatures(i); % 使用真实温度
    
    % 计算大气中CO2的分压
    P_CO2 = CO2_air_ppm * P_total / 1e6;
    
    % 计算CO2在水中的浓度
    CO2_sea = k_H * P_CO2;
    
    % 定义微分方程，现在也接受温度T作为参数
    odefun_with_CO2 = @(t, y) npzd_ode_with_CO2(t, y, Um, kN, Gm, lambda, gamma, theta, Is, I0, Qg10, Qh10, Mp, MZ, e, CO2_sea, k_CO2, CO2_max, CO2_air_ppm, alpha, P_CO2, T);
    
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
    if max(P_last_section) <= 0.3 * P_initial
        crash_CO2 = [crash_CO2, CO2_air_ppm]; % 记录崩溃时的CO2浓度
    else
        no_crash_CO2 = [no_crash_CO2, CO2_air_ppm]; % 记录未崩溃的CO2浓度
    end
end

% 排序未崩溃的CO2浓度（按从高到低排序）
sorted_no_crash_CO2 = sort(no_crash_CO2, 'descend'); % 按从高到低排序

% 输出没有崩溃的CO2浓度及对应的温度值
fprintf('CO2浓度 (ppm)  对应温度 (°C)\n');
fprintf('------------------------------\n');
for i = 1:length(sorted_no_crash_CO2)
    idx = find(realCO2Concentrations == sorted_no_crash_CO2(i));
    corresponding_temp = realTemperatures(idx);
    fprintf('%9.2f        %9.2f\n', sorted_no_crash_CO2(i), corresponding_temp);
end


% 定义微分方程
function dydt = npzd_ode_with_CO2(t, y, Um, kN, Gm, lambda, gamma, theta, Is, I0, Qg10, Qh10, Mp, MZ, e, CO2_concentration, k_CO2, C_max, CO2_air_ppm, alpha, acidification, T)
    N = y(1);
    P = y(2);
    Z = y(3);
    D = y(4);

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