clear; 
clc;

% 读取数据
data = readtable('fault_with_params1117-1.xlsx');

% 固定参数
mu = 10000000000; % [Pa]
nu = 0.2;
R = 6371000; % [m]
ert = 0;
etz = 0;
f = 0.75;

% 固定参考点
lon0 = 105.33; % 固定参考经度
lat0 = 29.235; % 固定参考纬度

% 结果存储
coulomb_stress = zeros(height(data), 1);
sigma_normal = zeros(height(data), 1);
tau_shear = zeros(height(data), 1);

% 遍历每一行数据
for i = 1:height(data)
    % 获取输入参数
    lonr = data{i, 2};
    latr = data{i, 1};
    ezz = data{i, 8};
    err = data{i, 9};
    ett = data{i, 10};
    ezr = data{i, 11};
    p = data{i, 12};
    st = data{i, 5};
    di = data{i, 6};
    ra = data{i, 7};
    
    % 计算两点的相对坐标x和y
    latf = deg2rad(latr);
    lonf = deg2rad(lonr);
    latt = deg2rad(lat0);
    lont = deg2rad(lon0);
    
    b = 0.5 * pi - latf;
    c = 0.5 * pi - latt;
    if lont > lonf
        aa = lont - lonf;
        if aa <= pi
            iangle = 1;
        else
            aa = 2 * pi - aa;
            iangle = -1;
        end
    else
        aa = lonf - lont;
        if aa <= pi
            iangle = -1;
        else
            aa = 2 * pi - aa;
            iangle = 1;
        end
    end
    s = cos(b) * cos(c) + sin(b) * sin(c) * cos(aa);
    a = acos(sign(min(abs(s), 1), s));
    dis = a * R;
    if a * b * c == 0
        angleb = 0;
        anglec = 0;
    else
        s = 0.5 * (a + b + c);
        anglec = 2 * asin(min(1, sqrt(sin(s - a) * sin(s - b) / (sin(a) * sin(b)))));
        angleb = 2 * asin(min(1, sqrt(sin(s - a) * sin(s - c) / (sin(a) * sin(c)))));
        if iangle == 1
            angleb = 2 * pi - angleb;
        else
            anglec = 2 * pi - anglec;
        end
    end
    x = dis * cos(anglec);
    y = dis * sin(anglec);
    azi = atan2d(y, x) + 180.0;

    % 计算转换到直角坐标系下的应力分量
    sigma_rr = 2 * mu / (1 - 2 * nu) * ((1 - nu) * err + nu * (ett + ezz));
    sigma_tt = 2 * mu / (1 - 2 * nu) * ((1 - nu) * ett + nu * (err + ezz));
    sigma_zz = 2 * mu / (1 - 2 * nu) * ((1 - nu) * ezz + nu * (ett + err));
    sigma_zr = 2 * mu * ezr;
    sigma_rt = 2 * mu * ert;
    sigma_tz = 2 * mu * etz;
    sigma = [sigma_rr, sigma_rt, sigma_zr;
             sigma_rt, sigma_tt, sigma_tz;
             sigma_zr, sigma_tz, sigma_zz];
    cs = cosd(azi);
    ss = sind(azi);
    exyz = [cs, ss, 0;
            -ss, cs, 0;
            0, 0, 1];
    swap = zeros(3, 3);
    scar = zeros(3, 3);
    for m = 1:3
        for n = 1:3
            swap(m, n) = 0.0;
            for k = 1:3
                swap(m, n) = swap(m, n) + sigma(m, k) * exyz(k, n);
            end
        end
    end
    for m = 1:3
        for n = 1:3
            scar(m, n) = 0.0;
            for k = 1:3
                scar(m, n) = scar(m, n) + exyz(m, k) * swap(k, n);
            end
        end
    end
    snn = scar(1, 1);
    see = scar(2, 2);
    sdd = scar(3, 3);
    sne = scar(1, 2);
    sed = scar(2, 3);
    sdn = scar(3, 1);
    
    % 计算库伦应力
    st0 = deg2rad(st);
    di0 = deg2rad(di);
    ra0 = deg2rad(ra);
    s = [snn, sne, sdn;
         sne, see, sed;
         sdn, sed, sdd];
    ns = [sin(di0) * cos(st0 + 0.5 * pi);
          sin(di0) * sin(st0 + 0.5 * pi);
          -cos(di0)];
    rst = [cos(st0);
           sin(st0);
           0.0];
    rdi = [cos(di0) * cos(st0 + 0.5 * pi);
           cos(di0) * sin(st0 + 0.5 * pi);
           sin(di0)];
    ts = rst * cos(ra0) - rdi * sin(ra0);
    sig = 0.0;
    tau = 0.0;
    for m = 1:3
        for n = 1:3
            sig = sig + ns(m) * s(m, n) * ns(n);
            tau = tau + ts(m) * s(m, n) * ns(n);
        end
    end
    cfs = tau + f * (sig + p);
    
    % 保存结果
    coulomb_stress(i) = cfs;
    sigma_normal(i) = sig;
    tau_shear(i) = tau;
end

% 将结果添加到表格中
data.Coulomb_Stress = coulomb_stress;
data.Sigma_Normal = sigma_normal;
data.Tau_Shear = tau_shear;

% 保存结果到文件
writetable(data, 'fault_with_params_result1117-1.xlsx');

% 辅助函数
function sgn = sign(a, b)
    if b >= 0
        sgn = a;
    else
        sgn = -a;
    end
end
