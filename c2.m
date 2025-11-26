%%本代码为范晓冉2023年写，2024年整理；功能为将'stress_results_gpu.csv'文件中  %%
%%柱坐标系下坐标z,d以及应力分量sigma_zz, sigma_rr, sigma_theta_theta, sigma_zr%%
%%插值到直角坐标系下（仅转坐标）。得到的结果为：x,y,z，'theta', 'a', 'b', 'c', 'e'，保存在%%
%%'final_output_with_theta.csv'文件中%%

clear;
% 读取数据，跳过第一行
opts = detectImportOptions('stress_results1101.csv');
opts.DataLines = [2 Inf];  % 从第二行开始读取
data = readtable('stress_results1101.csv', opts);
data.Properties.VariableNames = {'z', 'd', 'a', 'b', 'c', 'e', 'p'};

% 初始化用于存储各个 z 值组的单元数组
groups = cell(51, 1);  % 从 0 到 8000，共 81 个 z 值

% 按照 z 值分组
for i = 0:12000
    if any(data.z == i)
        groups{i/240 + 1} = data(data.z == i, :);
    end
end

% 此时，groups 数组中的每个元素都包含相同 z 值的数据行
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 对每个组进行处理
% 对每个组进行处理
for i = 1:length(groups)
    % 获取当前组的数据
    group_data = groups{i};

    % 检查该组是否为空
    if isempty(group_data)
        continue;  % 如果为空，则跳过当前循环
    end

    % 获取当前组的 z 值
    z_value = group_data.z(1);

    % 初始化临时输出数组
    temp_output = [];

    % 遍历 x 和 y，步长为 100
    for x = -12000:240:12000
        for y = -12000:240:12000
            % 计算 d
            d = sqrt(x^2 + y^2);

            % 插值 a, b
            a = interp1(group_data.d, group_data.a, d, 'linear', 'extrap');
            b = interp1(group_data.d, group_data.b, d, 'linear', 'extrap');
            c = interp1(group_data.d, group_data.c, d, 'linear', 'extrap');
            e = interp1(group_data.d, group_data.e, d, 'linear', 'extrap');
            p = interp1(group_data.d, group_data.p, d, 'linear', 'extrap');
            % 计算极角 theta
            theta = atan2(y, x)-pi/2;
            if theta <0
                theta = theta+2*pi;
            end

            % 将结果添加到临时输出数组
            temp_output = [temp_output; z_value, x, y, theta, a, b, c, e, p];
        end
    end

    % 将临时输出结果合并到最终输出
    if i == 1
        final_output = temp_output;
    else
        final_output = [final_output; temp_output];
    end
end

% 输出最终结果到文件
writetable(array2table(final_output, 'VariableNames', {'z', 'x', 'y', 'theta', 'a', 'b', 'c', 'e','p'}), 'final_output_with_theta1101.csv');

