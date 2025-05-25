N = 1000; %��Ԫ��
U = zeros(3,N);
F = zeros(3,N);
gamma = 2;
CFL = 0.8;
dx = 1.0 / N; %-0.5~0.5
tmax = 0.2;

t = 0;
x = linspace(-0.5,0.5,N);

% ѡ��TVD������������FVS����
TVDLimitersfun = @Minmod_limiter;
FVSfun = @FluxStegerWarming;

Lfunc = @(U) L(U, gamma, dx, FVSfun, TVDLimitersfun);

k = floor(N / 2);
U(:, 1:k) = repmat([1; 0; 1], 1, k);
U(:, k+1:N) = repmat([0.125; 0; 0.1], 1, N - k);

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
    
    % ʵʱ����ͼ��
    plot_result(U, h, gamma);
    drawnow;
end
