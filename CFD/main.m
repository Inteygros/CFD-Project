N = 1000; %��Ԫ��
U = zeros(3,N);
F = zeros(3,N);
gamma = 1.4;
CFL = 0.9;
dx = 1.0 / N; %-0.5~0.5
tmax = 0.2;

t = 0;
x = linspace(-0.5,0.5,N);

% ѡ�񷽷�
Option1 = 1; % 1��ʹ��TVD
                     % 2��ʹ��Ⱥ�ٶȿ���
                     % 3��ʹ��WENO
Option2 = 1; % 1��ʹ��FVS
                     % 2��ʹ��FDS
show_exact_solution = false; %�Ƿ���ʾ��ȷ��
                     
TVDLimitersfun = @Minmod_limiter; % TVD������������Minmod_limiter��/��MC_limiter��/��Superbee_limiter��
GVCscheme = @CD2_UD2; % GVC��ʽ���������Ĳ��+����ӭ�磨CD2_UD2��
WENOscheme = @WENO_5; % WENO��ʽ��5��WENO��WENO_5��
FVSfun = @Flux_StegerWarming; % FVS������Steger-Warming��Flux_StegerWarming��
UsingCharacteristicReconstruction = true; % FVS�Ƿ�ʹ�������ع�
FDSfun = @Flux_Roe; % FDS������Roe��ʽ��Flux_Roe��

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
    
    % ����RKʱ���ƽ�
    U = rk3_step(U, dt, Lfunc);
    U = real(U);
    Uexact = sod_exact_solution(t, N, gamma);
    % ʵʱ����ͼ��
    plot_result(U, h, gamma, t, N, show_exact_solution);
    drawnow;
end
