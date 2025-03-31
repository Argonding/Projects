% 定义变量
delta_C_w = -5;
delta_C_TOC = -25;
delta_C_carb = -12;
Delta = 25;

% 遍历 f_org
f_org_values = linspace(0, 1, 50);
f_ac_values = zeros(size(f_org_values));

for i = 1:length(f_org_values)
    f_org = f_org_values(i);
    
    % 计算 delta_C_prec 和 delta_C_org
    delta_C_prec = sym('x');
    delta_C_org = delta_C_prec - Delta;
    
    % 定义符号变量
    f_ac = sym('f_ac');
    delta_C_ac = sym('delta_C_ac');
    
    % 方程组
    eq1 = -5 == (1 - f_ac - f_org)*delta_C_prec + f_ac*delta_C_ac + f_org*(delta_C_prec - Delta);
    eq2 = -25 == (f_ac/(f_ac + f_org))*delta_C_ac + (f_org/(f_ac + f_org))*(delta_C_prec - Delta);
    eq3 = -12 == (f_ac/(1 - f_org))*delta_C_ac + ((1 - f_org - f_ac)/(1 - f_org))*delta_C_prec;
    
    % 求解方程组
    sol = solve([eq1, eq2, eq3], [f_ac, delta_C_ac, delta_C_prec]);
    
    if ~isempty(sol.f_ac)
        f_ac_values(i) = double(sol.f_ac);
    else
        f_ac_values(i) = NaN;
    end
end

% 绘制 fac-forg 图
figure;
plot(f_org_values, f_ac_values, 'LineWidth', 2);
xlabel('f_{org}');
ylabel('f_{ac}');
title('f_{ac} vs f_{org}');
grid on;