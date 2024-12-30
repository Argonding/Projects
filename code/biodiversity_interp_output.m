% 读取Excel文件
filename = 'pbdb_biodiversity_data_from_540Ma.xlsx';
data = readtable(filename);

% 提取所需的列
max_ma = data.max_ma;        % 最大年龄（开始时间）
min_ma = data.min_ma;        % 最小年龄（结束时间）
sampled_in_bin = data.sampled_in_bin; % 物种多样性数据

% 创建一个新的时间轴，用于绘图
geological_time = (max_ma + min_ma) / 2; % 以每个区间的中点作为地质时间

% 只保留420 Ma到0 Ma之间的数据
filter = geological_time >= 0 & geological_time <= 420;
geological_time_filtered = geological_time(filter);
sampled_in_bin_filtered = sampled_in_bin(filter);

% 创建新的时间点，每0.5 Ma进行插值
new_time_points = 0:0.5:420; % 从0到420 Ma，每0.5 Ma一个点（正序）

% 使用插值函数对原始数据进行插值
smoothed_biodiversity_interp = interp1(geological_time_filtered, sampled_in_bin_filtered, new_time_points, 'linear', 'extrap');

% 将时间和生物多样性数据合并成一个表格
output_table = table(new_time_points', smoothed_biodiversity_interp', 'VariableNames', {'GeologicalTime_Ma', 'SmoothedBiodiversity'});

% 将数据写入新的Excel文件
output_filename = 'dataset_CO2_biodiv.xlsx';
writetable(output_table, output_filename);

% 显示输出的路径
disp(['Data saved to: ', output_filename]);
