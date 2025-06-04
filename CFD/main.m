N = 1000; %��Ԫ��
U = zeros(3,N);
F = zeros(3,N);
gamma = 1.4;
CFL = 0.8;
dx = 1.0 / N; %-0.5~0.5
tmax = 0.2;

t = 0;
x = linspace(-0.5,0.5,N);

% ѡ�񷽷�
Option1 = 1; % 1��ʹ��TVD
                     % 2��ʹ��Ⱥ�ٶȿ���
                     % 3��ʹ��WENO
Option2 = 2; % 1��ʹ��FVS
                     % 2��ʹ��FDS
                     
TVDLimitersfun = @Minmod_limiter; % TVD������������Minmod_limiter��/��MC_limiter��/��Superbee_limiter��
GVCscheme = @CD2_UD2; % GVC��ʽ���������Ĳ��+����ӭ�磨CD2_UD2�� / ��Ͻ��¸�ʽ��UC3_SC4��
WENOscheme = @WENO_5; % WENO��ʽ�������뵱ǰֻ�䱸��5��WENO��WENO_5��
FVSfun = @Flux_StegerWarming; % FVS������������ֻ�䱸��Steger-Warming��Flux_StegerWarming��
UsingCharacteristicReconstruction = true; % FVS�Ƿ�ʹ�������ع�
FDSfun = @Flux_Roe; % FDS������������ֻ�䱸��Roe��Flux_Roe��%�ƻ��Ľ�һ��Rusanov


Lfunc = @(U) LU(U, gamma, dx, FVSfun, FDSfun, TVDLimitersfun, GVCscheme, WENOscheme, UsingCharacteristicReconstruction, Option1, Option2);

k = floor(N / 2);
U(:, 1:k) = repmat([1; 0; 1/(gamma-1)], 1, k);
U(:, k+1:N) = repmat([0.125; 0; 0.1/(gamma-1)], 1, N - k);

[rho, u, p, ~] = conservative_to_primitive(U, gamma);
hFig = figure('Position', [400, 200, 600, 400]);
h = init_plot(x, rho, u, p);
drawnow;

while(t<tmax)
    [rho, u, p, c] = conservative_to_primitive(U, gamma);
    dt = CFL * dx / max(abs(u) + c);
    t = t+dt;
    
    % ����RKʱ���ƽ�
    U = rk3_step(U, dt, Lfunc);
    U = real(U);
    
    % ʵʱ����ͼ��
    plot_result(U, h, gamma);
    drawnow;
end
