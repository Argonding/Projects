% Define the file name
filename = 'RSE.xlsx';

% Read the data from the Excel file, specifically rows 3 to 39 and columns A to K
data = readcell(filename, 'Range', 'A3:K39');

% Extract the 'position' column (which is the 2nd column in the 'data' variable now)
position = cell2mat(data(:, 2));

% Extract the data for the desired columns (V through Cu-EF)
data_numeric = cell2mat(data(:, 3:11));

V = data_numeric(:, 1);
Mo = data_numeric(:, 2);
Ni = data_numeric(:, 3);
Cu = data_numeric(:, 4);
Al2O3 = data_numeric(:, 5);
V_EF = data_numeric(:, 6);
Mo_EF = data_numeric(:, 7);
Ni_EF = data_numeric(:, 8);
Cu_EF = data_numeric(:, 9);

% Combine data and position for sorting
combined_data = [position, V, Mo, Ni, Cu, Al2O3, V_EF, Mo_EF, Ni_EF, Cu_EF];

% Sort the data by 'position' in descending order
sorted_combined_data = sortrows(combined_data, 1, 'descend');

% Separate the sorted data
sorted_position = sorted_combined_data(:, 1);
sorted_V = sorted_combined_data(:, 2);
sorted_Mo = sorted_combined_data(:, 3);
sorted_Ni = sorted_combined_data(:, 4);
sorted_Cu = sorted_combined_data(:, 5);
sorted_Al2O3 = sorted_combined_data(:, 6);
sorted_V_EF = sorted_combined_data(:, 7);
sorted_Mo_EF = sorted_combined_data(:, 8);
sorted_Ni_EF = sorted_combined_data(:, 9);
sorted_Cu_EF = sorted_combined_data(:, 10);


% Define the custom color palette (RGB values normalized to [0, 1])
customColors = [
    66  52  88;    % Color 1
    119 96  142;   % Color 2
    161 142 174;   % Color 3
    35  51  66;    % Color 4
    74  90  105;   % Color 5
    167 174 190;   % Color 6
    84  60  48;    % Color 7
    151 132 115;   % Color 8
    184 150 122    % Color 9
] / 255; % Normalize RGB values by dividing by 255

% Create a figure and set its position and size
figure('Position', [100, 100, 1600, 300]); % Wide and relatively short figure

% Define common plot properties for consistency
plot_line_spec = '-o'; % Lines with circle markers
marker_size = 4;      % Smaller marker size
font_size_title = 10;
font_size_label = 9;

% Get the overall Y-axis limits from the sorted position data
overall_y_min_data = min(sorted_position);
overall_y_max_data = max(sorted_position);

% --- Create a background axes that spans the entire figure ---
% This axes will be invisible and used only for drawing the background lines
h_bg_axes = axes('Position', [0 0 1 1], ...         % [left, bottom, width, height] = spans entire figure
                 'XLim', [0 1], ...                 % Normalized X-limits (0 to 1 across the figure)
                 'YLim', [overall_y_min_data, overall_y_max_data], ... % Y-limits match the data range
                 'Visible', 'off', ...              % Make this background axes invisible
                 'Tag', 'BackgroundLines');         % Give it a tag for easier identification (optional)

% Send this background axes to the bottom layer so subplots draw on top
uistack(h_bg_axes, 'bottom');

% Draw the first horizontal line on the background axes
plot(h_bg_axes, [0 1], [line_y1, line_y1], ... % X-coordinates 0 to 1 (full figure width)
     'Color', line_color, ...
     'LineStyle', line_style, ...
     'LineWidth', line_width);

% Draw the second horizontal line on the background axes
plot(h_bg_axes, [0 1], [line_y2, line_y2], ... % X-coordinates 0 to 1 (full figure width)
     'Color', line_color, ...
     'LineStyle', line_style, ...
     'LineWidth', line_width);

% --- Now, create and plot the subplots as before (these will draw on top) ---
% Data for plotting (for easier iteration and consistent labeling)
plot_data_x = {sorted_V, sorted_Mo, sorted_Ni, sorted_Cu, sorted_Al2O3, ...
               sorted_V_EF, sorted_Mo_EF, sorted_Ni_EF, sorted_Cu_EF};
plot_titles = {'V', 'Mo', 'Ni', 'Cu', 'Al_{2}O_{3}', 'V-EF', 'Mo-EF', 'Ni-EF', 'Cu-EF'}; % Al2O3 title changed to LaTeX for subscript
plot_xlabels = {'Conc.', 'Conc.', 'Conc.', 'Conc.', 'Conc.', 'Factor', 'Factor', 'Factor', 'Factor'};

% 设置阴影区域的范围
shade_y1 = 2745;
shade_y2 = 3100;
shade_color = [0.6 0.6 0.6]; % 浅灰色

% Loop through all 9 subplots
for i = 1:9
    subplot(1,9,i); % 激活当前子图
    
    % 绘制数据线（先画，才能让xlim自动调整）
    plot(plot_data_x{i}, sorted_position, plot_line_spec, ...
         'Color', customColors(i,:), ...
         'MarkerFaceColor', customColors(i,:), ...
         'MarkerSize', marker_size);
    hold on;

    % 获取当前子图的 x 轴范围（这时已由 plot 自动设定好）
    x_limits = xlim();

    % 绘制阴影区域，覆盖整个横轴
    fill([x_limits(1), x_limits(2), x_limits(2), x_limits(1)], ...
         [shade_y1, shade_y1, shade_y2, shade_y2], ...
         shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % 再次绘制数据线，确保在阴影上层
    plot(plot_data_x{i}, sorted_position, plot_line_spec, ...
         'Color', customColors(i,:), ...
         'MarkerFaceColor', customColors(i,:), ...
         'MarkerSize', marker_size);

    % 设置标题和标签
    title(plot_titles{i}, 'FontSize', font_size_title)

    if i == 1
        ylabel('Position', 'FontSize', font_size_label);
    end

    set(gca, 'YDir', 'normal'); 
    grid on;
end



% Add a super title for the entire figure
sgtitle('RSE concentrations around F-F boundary');