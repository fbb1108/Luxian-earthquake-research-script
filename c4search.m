clear;

% 读取final_output_with_lat_lon.csv
final_data = readtable('final_output_with_lat_lon1101.csv');

% 读取fault.xlsx
fault_data = readtable('fault2-2.xlsx');

% 获取final_output_with_lat_lon.csv的列
final_depth = final_data{:, 1};
final_lat = final_data{:, 10};
final_lon = final_data{:, 11};
final_params = final_data{:, 5:9};

% 获取fault.xlsx的列
fault_depth = fault_data{:, 3};
fault_lat = fault_data{:, 1};
fault_lon = fault_data{:, 2};

% 初始化保存结果的数组
num_fault_rows = size(fault_data, 1);
new_params = zeros(num_fault_rows, 5);
closest_lat = zeros(num_fault_rows, 1);
closest_lon = zeros(num_fault_rows, 1);

% 搜索符合条件的行
for i = 1:num_fault_rows
    % 搜索最近的经度和纬度
    distances = sqrt((final_lat - fault_lat(i)).^2 + (final_lon - fault_lon(i)).^2);
    [~, min_idx] = min(distances);
    
    % 获取最接近经纬度的所有行
    matched_lat = final_lat(min_idx);
    matched_lon = final_lon(min_idx);
    matched_rows = final_data(final_lat == matched_lat & final_lon == matched_lon, :);
    
    % 记录最接近的经纬度
    closest_lat(i) = matched_lat;
    closest_lon(i) = matched_lon;
    
    % 找到深度上最接近的两个点
    [sorted_depths, sort_idx] = sort(matched_rows{:, 1});
    matched_rows = matched_rows(sort_idx, :);
    
    if length(sorted_depths) < 2
        % 如果没有足够的点进行插值，跳过该行
        new_params(i, :) = NaN;
        continue;
    end
    
    % 找到深度最接近的两个点
    lower_idx = find(sorted_depths <= fault_depth(i), 1, 'last');
    upper_idx = find(sorted_depths >= fault_depth(i), 1, 'first');
    
    if isempty(lower_idx) || isempty(upper_idx)
        new_params(i, :) = NaN;
        continue;
    end
    
    if lower_idx == upper_idx
        % 如果深度正好匹配
        new_params(i, :) = matched_rows{lower_idx, 5:9};
    else
        % 进行线性插值
        lower_depth = sorted_depths(lower_idx);
        upper_depth = sorted_depths(upper_idx);
        
        lower_params = matched_rows{lower_idx, 5:9};
        upper_params = matched_rows{upper_idx, 5:9};
        
        ratio = (fault_depth(i) - lower_depth) / (upper_depth - lower_depth);
        new_params(i, :) = lower_params + ratio * (upper_params - lower_params);
    end
end

% 将新的参数和最接近的经纬度添加到fault_data中
fault_data = [fault_data, array2table(new_params, 'VariableNames', {'Param5', 'Param6', 'Param7', 'Param8', 'Param9'}), ...
              array2table(closest_lat, 'VariableNames', {'ClosestLat'}), ...
              array2table(closest_lon, 'VariableNames', {'ClosestLon'})];

% 保存结果到新的Excel文件
writetable(fault_data, 'fault_with_params1117-2.xlsx');
