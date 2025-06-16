N = 1000; %单元数
U = zeros(3,N);
F = zeros(3,N);
gamma = 1.4;
CFL = 0.9;
dx = 1.0 / N; %-0.5~0.5
tmax = 0.2;

t = 0;
x = linspace(-0.5,0.5,N);

% 选择方法
Option1 = 1; % 1：使用TVD
                     % 2：使用群速度控制
                     % 3：使用WENO
Option2 = 1; % 1：使用FVS
                     % 2：使用FDS
show_exact_solution = false; %是否显示精确解
                     
TVDLimitersfun = @Minmod_limiter; % TVD限制器函数（Minmod_limiter）/（MC_limiter）/（Superbee_limiter）
GVCscheme = @CD2_UD2; % GVC格式：二阶中心差分+二阶迎风（CD2_UD2）
WENOscheme = @WENO_5; % WENO格式：5阶WENO（WENO_5）
FVSfun = @Flux_StegerWarming; % FVS函数：Steger-Warming（Flux_StegerWarming）
UsingCharacteristicReconstruction = true; % FVS是否使用特征重构
FDSfun = @Flux_Roe; % FDS函数：Roe格式（Flux_Roe）

k = floor(N / 2);
U(:, 1:k) = repmat([1; 0; 1/(gamma-1)], 1, k);
U(:, k+1:N) = repmat([0.125; 0; 0.1/(gamma-1)], 1, N - k);
[~, u, ~, c] = conservative_to_primitive(U, gamma);
dt = CFL * dx / max(abs(u) + c);

Lfunc = @(U, dt) LU(U, gamma, dx, dt, FVSfun, FDSfun, TVDLimitersfun, GVCscheme, WENOscheme, UsingCharacteristicReconstruction, Option1, Option2);
    
hFig = figure('Position', [400, 200, 600, 400]);
h = init_plot(x, U, gamma, show_exact_solution);
drawnow;

while(t<tmax)
    [~, u, ~, c] = conservative_to_primitive(U, gamma);
    dt = CFL * dx / max(abs(u) + c);
    t = t+dt;
    
    % 三阶RK时间推进
    U = rk3_step(U, dt, Lfunc);
    U = real(U);
    Uexact = sod_exact_solution(t, N, gamma);
    % 实时更新图像
    plot_result(U, h, gamma, t, N, show_exact_solution);
    drawnow;
end
