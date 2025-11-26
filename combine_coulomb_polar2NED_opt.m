function combine_coulomb_polar2NED_opt(file1, ref1, file2, ref2, outfile)
% 线性叠加两套以“柱坐标(r–t–z)”给出的应变/孔压场对同一断层的影响
% 额外新增：分别对 file1 场、file2 场、叠加场，搜索最优断层面参数（含两组共轭解）
%
% 输入：
%   file1, file2 : 两个 xlsx（行列一致，前7列一致）
%   ref1 = [lon0_1, lat0_1] ：文件1使用的参考点
%   ref2 = [lon0_2, lat0_2] ：文件2使用的参考点
%   outfile : 输出文件名
%
% 列定义（按你的代码）：
%   1: lat, 2: lon, 3: depth_m, 4: depth_km, 5: strike, 6: dip, 7: rake,
%   8: ezz, 9: err, 10: ett, 11: ezr, 12: p
% 注意：没有第13列的 ert；本脚本中固定 ε_NE=0、ε_ED=0、ε_tz=0

mu = 10e9;   % Pa
nu = 0.2;
f  = 0.4;    % 摩擦系数

% 搜索参数（可按需调整）
dStrikeCoarse = 5;    % 粗扫 走向步长（度）
dDipCoarse    = 2;    % 粗扫 倾角步长（度）
refineWinS    = 6;    % 精化窗口 ±度数（走向）
refineWinD    = 4;    % 精化窗口 ±度数（倾角）
dStrikeFine   = 0.5;  % 精扫 走向步长（度）
dDipFine      = 0.5;  % 精扫 倾角步长（度）
minSepDeg     = 25;   % 两个最优面法向间的最小分离角（度），避免重复近邻解

% 列索引（按你的文件）
COL.lat = 1; COL.lon = 2; COL.strike = 5; COL.dip = 6; COL.rake = 7;
COL.ezz = 8; COL.err = 9; COL.ett = 10; COL.ezr = 11; COL.p = 12;

T1 = readtable(file1);
T2 = readtable(file2);
assert(height(T1)==height(T2), '两文件行数不同');
assert(all(T1{:,1:7}==T2{:,1:7}, 'all'), '前7列(位置/几何)不一致');

n = height(T1);

% 预分配输出
sn_tot = zeros(n,1);  ta_tot = zeros(n,1);  cfs_tot = zeros(n,1);
sn1    = zeros(n,1);  ta1    = zeros(n,1);  cfs1    = zeros(n,1);
sn2    = zeros(n,1);  ta2    = zeros(n,1);  cfs2    = zeros(n,1);
p_out1 = zeros(n,1);  p_out2 = zeros(n,1);  p_outT  = zeros(n,1);

% ——最优面输出（每场两组解）——
opt1_f1 = init_opt_arrays(n);   % file1 场
opt2_f1 = init_opt_arrays(n);

opt1_f2 = init_opt_arrays(n);   % file2 场
opt2_f2 = init_opt_arrays(n);

opt1_ft = init_opt_arrays(n);   % 叠加场
opt2_ft = init_opt_arrays(n);

for i = 1:n
    st = T1{i, COL.strike}; di = T1{i, COL.dip}; ra = T1{i, COL.rake};
    [nvec_given, tvec_given] = fault_basis_NED(st, di, ra);

    % --- 文件1 ---
    az1      = local_azimuth_from_ref(T1{i,COL.lon}, T1{i,COL.lat}, ref1(1), ref1(2));
    S1_NED   = build_stress_in_NED_from_polar(T1(i,:), COL, az1, mu, nu);
    p1       = T1{i, COL.p};

    % --- 文件2 ---
    az2      = local_azimuth_from_ref(T2{i,COL.lon}, T2{i,COL.lat}, ref2(1), ref2(2));
    S2_NED   = build_stress_in_NED_from_polar(T2(i,:), COL, az2, mu, nu);
    p2       = T2{i, COL.p};

    % ---- 线性叠加 ----
    S_tot = S1_NED + S2_NED;
    p_tot = p1 + p2;

    % ---- 给定断层上的单场 / 合成 CFS（与你原逻辑一致）----
    sn1(i) = nvec_given' * S1_NED * nvec_given;
    ta1(i) = tvec_given' * S1_NED * nvec_given;
    cfs1(i) = ta1(i) + f * (sn1(i) + p1);

    sn2(i) = nvec_given' * S2_NED * nvec_given;
    ta2(i) = tvec_given' * S2_NED * nvec_given;
    cfs2(i) = ta2(i) + f * (sn2(i) + p2);

    sn_tot(i)  = nvec_given' * S_tot * nvec_given;
    ta_tot(i)  = tvec_given' * S_tot * nvec_given;
    cfs_tot(i) = ta_tot(i) + f * (sn_tot(i) + p_tot);

    % ---- 记录孔压 ----
    p_out1(i) = p1;  p_out2(i) = p2;  p_outT(i) = p_tot;

    % ===== 搜索最优断层面（file1、file2、叠加场） =====
    optPair_f1 = search_optimal_two(S1_NED, p1, f, dStrikeCoarse, dDipCoarse, ...
                                    refineWinS, refineWinD, dStrikeFine, dDipFine, minSepDeg);
    optPair_f2 = search_optimal_two(S2_NED, p2, f, dStrikeCoarse, dDipCoarse, ...
                                    refineWinS, refineWinD, dStrikeFine, dDipFine, minSepDeg);
    optPair_ft = search_optimal_two(S_tot , p_tot, f, dStrikeCoarse, dDipCoarse, ...
                                    refineWinS, refineWinD, dStrikeFine, dDipFine, minSepDeg);

    opt1_f1 = assign_opt_row(opt1_f1, i, optPair_f1.best);
    opt2_f1 = assign_opt_row(opt2_f1, i, optPair_f1.second);

    opt1_f2 = assign_opt_row(opt1_f2, i, optPair_f2.best);
    opt2_f2 = assign_opt_row(opt2_f2, i, optPair_f2.second);

    opt1_ft = assign_opt_row(opt1_ft, i, optPair_ft.best);
    opt2_ft = assign_opt_row(opt2_ft, i, optPair_ft.second);
end

% 写出（保留原列 + 新列）
Tout = T1;
Tout.pore_pressure_1    = p_out1;   % 第一组孔压 (Pa)
Tout.pore_pressure_2    = p_out2;   % 第二组孔压 (Pa)
Tout.pore_pressure_tot  = p_outT;   % 合成孔压 (Pa)

Tout.tau_shear_1        = ta1;
Tout.sigma_normal_1     = sn1;
Tout.Coulomb_Stress_1   = cfs1;

Tout.tau_shear_2        = ta2;
Tout.sigma_normal_2     = sn2;
Tout.Coulomb_Stress_2   = cfs2;

Tout.tau_shear_tot      = ta_tot;
Tout.sigma_normal_tot   = sn_tot;
Tout.Coulomb_Stress_tot = cfs_tot;

% ——新增：最优面参数列（file1）——
Tout.stk_opt1_f1 = opt1_f1.stk;  Tout.dip_opt1_f1 = opt1_f1.dip;  Tout.rak_opt1_f1 = opt1_f1.rak;
Tout.cfs_opt1_f1 = opt1_f1.cfs;  Tout.sn_opt1_f1  = opt1_f1.sn;   Tout.tau_opt1_f1  = opt1_f1.tau;

Tout.stk_opt2_f1 = opt2_f1.stk;  Tout.dip_opt2_f1 = opt2_f1.dip;  Tout.rak_opt2_f1 = opt2_f1.rak;
Tout.cfs_opt2_f1 = opt2_f1.cfs;  Tout.sn_opt2_f1  = opt2_f1.sn;   Tout.tau_opt2_f1  = opt2_f1.tau;

% ——新增：最优面参数列（file2）——
Tout.stk_opt1_f2 = opt1_f2.stk;  Tout.dip_opt1_f2 = opt1_f2.dip;  Tout.rak_opt1_f2 = opt1_f2.rak;
Tout.cfs_opt1_f2 = opt1_f2.cfs;  Tout.sn_opt1_f2  = opt1_f2.sn;   Tout.tau_opt1_f2  = opt1_f2.tau;

Tout.stk_opt2_f2 = opt2_f2.stk;  Tout.dip_opt2_f2 = opt2_f2.dip;  Tout.rak_opt2_f2 = opt2_f2.rak;
Tout.cfs_opt2_f2 = opt2_f2.cfs;  Tout.sn_opt2_f2  = opt2_f2.sn;   Tout.tau_opt2_f2  = opt2_f2.tau;

% ——新增：最优面参数列（叠加场）——
Tout.stk_opt1_ft = opt1_ft.stk;  Tout.dip_opt1_ft = opt1_ft.dip;  Tout.rak_opt1_ft = opt1_ft.rak;
Tout.cfs_opt1_ft = opt1_ft.cfs;  Tout.sn_opt1_ft  = opt1_ft.sn;   Tout.tau_opt1_ft  = opt1_ft.tau;

Tout.stk_opt2_ft = opt2_ft.stk;  Tout.dip_opt2_ft = opt2_ft.dip;  Tout.rak_opt2_ft = opt2_ft.rak;
Tout.cfs_opt2_ft = opt2_ft.cfs;  Tout.sn_opt2_ft  = opt2_ft.sn;   Tout.tau_opt2_ft  = opt2_ft.tau;

writetable(Tout, outfile);
fprintf('Done. Wrote %s\n', outfile);
end

% ====== 工具：初始化最优结构数组 ======
function S = init_opt_arrays(n)
S.stk = zeros(n,1);
S.dip = zeros(n,1);
S.rak = zeros(n,1);
S.cfs = zeros(n,1);
S.sn  = zeros(n,1);
S.tau = zeros(n,1);
end

function S = assign_opt_row(S, i, O)
S.stk(i) = O.strike;
S.dip(i) = O.dip;
S.rak(i) = O.rake;
S.cfs(i) = O.cfs;
S.sn(i)  = O.sn;
S.tau(i) = O.tau;
end

% ====== 最优面两组解搜索（粗扫+精化） ======
function out = search_optimal_two(S, p, f, dS_coarse, dD_coarse, winS, winD, dS_fine, dD_fine, minSepDeg)
% 返回两个解：best / second（若找不到第二个，就复制best）
best   = struct('strike',NaN,'dip',NaN,'rake',NaN,'cfs',-Inf,'sn',NaN,'tau',NaN,'n',[NaN;NaN;NaN]);
second = best;

% --- 粗扫 ---
[best, pool] = grid_search_block(S, p, f, 0, 360-dS_coarse, dS_coarse, 1, 89, dD_coarse, best);

% 在 best 周边做精化
smin = max(0,   best.strike - winS);
smax = min(360, best.strike + winS);
dmin = max(1,   best.dip    - winD);
dmax = min(89,  best.dip    + winD);
best = grid_search_refine(S, p, f, smin, smax, dS_fine, dmin, dmax, dD_fine, best);

% 寻找第二个峰值：从 pool 中挑与 best 分离足够的候选，再对其做精化
cand = filter_by_separation(pool, best.n, minSepDeg);
if ~isempty(cand)
    [~, idx] = max([cand.cfs]);
    second0 = cand(idx);
    % 对 second0 周边精化
    smin = max(0,   second0.strike - winS);
    smax = min(360, second0.strike + winS);
    dmin = max(1,   second0.dip    - winD);
    dmax = min(89,  second0.dip    + winD);
    second = grid_search_refine(S, p, f, smin, smax, dS_fine, dmin, dmax, dD_fine, second0);
else
    second = best; % 若没有明显第二峰，则复制best
end

out.best = best;
out.second = second;
end

function [best, pool] = grid_search_block(S, p, f, s0, s1, ds, d0, d1, dd, best)
    % 初始化空结构体，用于保存所有候选面
    pool = struct('cfs', {}, 'sn', {}, 'tau', {}, 'strike', {}, 'dip', {}, 'rake', {}, 'n', {});
    for strike = s0:ds:s1
        for dip = d0:dd:d1
            [cfs, sn, tau, rake, nvec] = cfs_on_plane(S, p, f, strike, dip);
            if cfs > best.cfs
                best.cfs = cfs; best.sn = sn; best.tau = tau;
                best.strike = mod(strike,360); best.dip = dip; best.rake = rake; best.n = nvec;
            end
            pool(end+1) = struct('cfs', cfs, 'sn', sn, 'tau', tau, ...
                                 'strike', mod(strike,360), 'dip', dip, 'rake', rake, 'n', nvec);
        end
    end
end



function best = grid_search_refine(S, p, f, smin, smax, ds, dmin, dmax, dd, best)
for strike = smin:ds:smax
    for dip = dmin:dd:dmax
        [cfs, sn, tau, rake, nvec] = cfs_on_plane(S, p, f, strike, dip);
        if cfs > best.cfs
            best.cfs = cfs; best.sn = sn; best.tau = tau;
            best.strike = mod(strike,360); best.dip = dip; best.rake = rake; best.n = nvec;
        end
    end
end
end

function C = filter_by_separation(pool, n_best, minSepDeg)
    % 返回与 n_best 法向分离角 >= minSepDeg 的候选集合（结构体数组）
    if isempty(pool)
        C = struct('cfs', {}, 'sn', {}, 'tau', {}, 'strike', {}, 'dip', {}, 'rake', {}, 'n', {});
        return
    end
    C = struct('cfs', {}, 'sn', {}, 'tau', {}, 'strike', {}, 'dip', {}, 'rake', {}, 'n', {}); % 空结构体
    for k = 1:numel(pool)
        ang = vec_angle_deg(pool(k).n, n_best);
        if ang >= minSepDeg
            C(end+1) = pool(k); %#ok<AGROW>
        end
    end
end


function ang = vec_angle_deg(a, b)
a = a(:)/norm(a); b = b(:)/norm(b);
ang = rad2deg(acos( max(-1, min(1, dot(a,b))) ));
end

% ====== 给定面（strike,dip）上计算 CFS 与 rake（rake 取使剪应力最大的方向） ======
function [cfs, sig_n, tau, rake, nvec] = cfs_on_plane(S, p, f, strike, dip)
% 法向
nvec = normal_from_strike_dip(strike, dip);

% 面上牵引
t = S * nvec;
sig_n = dot(nvec, t);

% 面内剪应力向量与模
t_tan = t - sig_n * nvec;
tau = norm(t_tan);

% 本实现与主程序保持一致：ΔCFS = τ + f * (σ_n + p)
cfs = tau + f * (sig_n + p);

% 计算 rake：先构造走向/倾向基，再把 s = t_tan/||t_tan|| 展开
[sdir, ddir] = strike_dip_axes(strike, dip);
if tau < eps
    rake = 0; % 剪应力很小，rake 不敏感
else
    svec = t_tan / tau;
    % 与 fault_basis_NED 中关系保持：svec = sdir*cos(rake) - ddir*sin(rake)
    rake = atan2( -dot(svec, ddir), dot(svec, sdir) ); % 弧度
    rake = rad2deg(wrapToPi(rake));
end
end

function nvec = normal_from_strike_dip(strike_deg, dip_deg)
st = deg2rad(strike_deg);
di = deg2rad(dip_deg);
% 与你现有 fault_basis_NED 中 nvec 一致（Z向下）
nvec = [  sin(di)*cos(st+pi/2);
          sin(di)*sin(st+pi/2);
         -cos(di) ];
nvec = nvec / norm(nvec);
end

function [sdir, ddir] = strike_dip_axes(strike_deg, dip_deg)
st = deg2rad(strike_deg);
di = deg2rad(dip_deg);
sdir = [cos(st); sin(st); 0];                        % 走向方向（水平）
ddir = [cos(st+pi/2)*cos(di);  sin(st+pi/2)*cos(di);  sin(di)]; % 倾向方向（指向下倾）
end

% ====== 你已有的函数（保持不变） ======

function S_NED = build_stress_in_NED_from_polar(row, COL, az_deg, mu, nu)
ezz = row{1, COL.ezz};  err = row{1, COL.err};  ett = row{1, COL.ett};
ezr = row{1, COL.ezr};  ert = 0;                etz = 0;

sigma_rr = 2*mu/(1-2*nu)*((1-nu)*err + nu*(ett+ezz));
sigma_tt = 2*mu/(1-2*nu)*((1-nu)*ett + nu*(err+ezz));
sigma_zz = 2*mu/(1-2*nu)*((1-nu)*ezz + nu*(ett+err));
sigma_zr = 2*mu*ezr;    sigma_rt = 2*mu*ert;    sigma_tz = 2*mu*etz;

S_rtz = [sigma_rr, sigma_rt, sigma_zr;
         sigma_rt, sigma_tt, sigma_tz;
         sigma_zr, sigma_tz, sigma_zz];

az = deg2rad(az_deg); c = cos(az); s = sin(az);
R = [ c,  s,  0;   % (r,t,z) -> (N,E,D)
     -s,  c,  0;
      0,  0,  1];
S_NED = R * S_rtz * R.';
end

function az = local_azimuth_from_ref(lon, lat, lon0, lat0)
[az, ~] = azimuth_geodesic(lat0, lon0, lat, lon); % deg
end

function [az, dist] = azimuth_geodesic(lat1, lon1, lat2, lon2)
R = 6371000;
phi1 = deg2rad(lat1); lam1 = deg2rad(lon1);
phi2 = deg2rad(lat2); lam2 = deg2rad(lon2);
dlam = lam2 - lam1;
y = sin(dlam) * cos(phi2);
x = cos(phi1)*sin(phi2) - sin(phi1)*cos(phi2)*cos(dlam);
az = mod(rad2deg(atan2(y, x)), 360);
c = acos( max(-1,min(1, sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(dlam) )) );
dist = R * c;
end

function [nvec, tvec] = fault_basis_NED(strike_deg, dip_deg, rake_deg)
st = deg2rad(strike_deg); di = deg2rad(dip_deg); ra = deg2rad(rake_deg);
sdir = [cos(st); sin(st); 0];
ddir = [cos(st+pi/2)*cos(di);  sin(st+pi/2)*cos(di);  sin(di)];
nvec = [  sin(di)*cos(st+pi/2);
          sin(di)*sin(st+pi/2);
         -cos(di) ];
nvec = nvec / norm(nvec);
tvec = sdir*cos(ra) - ddir*sin(ra);
tvec = tvec / norm(tvec);
end
