% 读取数据
data = readmatrix('final_output_with_theta1101.csv');

% 提取深度 (z), x, 和 y 坐标
z = data(:, 1);
x = data(:, 2);
y = data(:, 3);

% 定义原点的纬度和经度 (以度为单位)
origin_lat = 29.235; % 替换为实际的原点纬度
origin_lon = 105.33; % 替换为实际的原点经度

% 地球半径 (单位: 米)
R = 6371000;

% 预分配纬度和经度数组
latitudes = zeros(size(x));
longitudes = zeros(size(x));

% 转换 x, y 到纬度和经度
for i = 1:length(x)
    % 水平距离转换为纬度
    dlat = rad2deg(y(i) / R);
    % 水平距离转换为经度
    dlon = rad2deg(x(i) / (R * cosd(origin_lat)));
    
    % 计算最终的纬度和经度
    latitudes(i) = origin_lat + dlat;
    longitudes(i) = origin_lon + dlon;
end

% 将新的纬度和经度附加到原始数据
new_data = [data, latitudes, longitudes];

% 将新数据写入一个新的CSV文件
writematrix(new_data, 'final_output_with_lat_lon1101.csv');
