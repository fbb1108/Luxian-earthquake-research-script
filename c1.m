clear;
% 读取数据文件，跳过第一行文字描述
dataFile = 'h79.dat';  
data = importdata(dataFile, ' ', 1);

% 提取数据列
z = data.data(:, 1);
d = data.data(:, 2);
ezz = data.data(:, 5);
err = data.data(:, 6);
ett = data.data(:, 7);
ezr = data.data(:, 8);
pp = data.data(:, 10);


% 创建保存结果的矩阵
stress_results = [z, d, ezz, err, ett, ezr, pp];

% 保存结果到文件
outputFile = 'stress_results1101.csv';  % 请替换为你想要保存结果的文件路径
csvwrite(outputFile, stress_results);
