% Sod问题精确解
function U = sod_exact_solution(t, N, gamma)
    x = linspace(-0.5, 0.5, N);
    
    % 初始状态
    rho_l = 1.0;
    u_l = 0.0;
    p_l = 1.0;
    rho_r = 0.125;
    u_r = 0.0;
    p_r = 0.1;
    
    % 声速计算
    c_l = sqrt(gamma * p_l / rho_l);
    c_r = sqrt(gamma * p_r / rho_r);
    
    % 求解黎曼问题的中间状态
    [p_star, u_star] = solve_riemann(rho_l, u_l, p_l, rho_r, u_r, p_r, gamma);
    
    % 中间状态密度
    rho_l_star = rho_l * (p_star/p_l)^(1/gamma); % 左侧稀疏波
    rho_r_star = rho_r * ( (gamma+1)*p_star/p_r + gamma - 1 ) / ...
                ( (gamma+1) + (gamma-1)*p_star/p_r ); % 右侧激波
    
    % 中间状态声速
    c_l_star = c_l * (p_star/p_l)^((gamma-1)/(2*gamma));
    
    % 波速计算
    s_head = u_l - c_l;                   % 稀疏波头
    s_tail = u_star - c_l_star;            % 稀疏波尾
    s_shock = u_r + c_r * sqrt(1 + (gamma+1)/(2*gamma)*(p_star/p_r - 1)); % 激波速度
    
    % 初始化变量
    rho = zeros(1, N);
    u = zeros(1, N);
    p = zeros(1, N);
    
    % 对于每个位置点计算精确解
    for i = 1:N
        xi = x(i) / t;  % 自相似坐标
        
        if xi <= s_head
            % 区域1: 左静止区
            rho(i) = rho_l;
            u(i) = u_l;
            p(i) = p_l;
            
        elseif xi <= s_tail
            % 区域2: 稀疏波区
            constant = u_l + 2*c_l/(gamma-1);
            u(i) = xi + (gamma-1)/(gamma+1)*(constant - xi);
            c = (gamma-1)/(gamma+1)*(constant - xi);
            rho(i) = rho_l * (c/c_l)^(2/(gamma-1));
            p(i) = p_l * (rho(i)/rho_l)^gamma;
            
        elseif xi <= u_star
            % 区域3: 接触间断左侧
            rho(i) = rho_l_star;
            u(i) = u_star;
            p(i) = p_star;
            
        elseif xi <= s_shock
            % 区域4: 接触间断右侧
            rho(i) = rho_r_star;
            u(i) = u_star;
            p(i) = p_star;
            
        else
            % 区域5: 右静止区
            rho(i) = rho_r;
            u(i) = u_r;
            p(i) = p_r;
        end
    end
    
    % 转换为守恒变量
    E = p/(gamma-1) + 0.5*rho.*u.^2;
    U = [rho; rho.*u; E];

    % 黎曼求解器
    function [p_star, u_star] = solve_riemann(rho_l, u_l, p_l, rho_r, u_r, p_r, gamma)
        % 牛顿迭代法
        max_iter = 20;
        p_guess = 0.3;
        
        for iter = 1:max_iter
            % 左侧稀疏波函数
            c_l = sqrt(gamma * p_l / rho_l);
            f_l = 2 * c_l / (gamma - 1) * ((p_guess/p_l)^((gamma-1)/(2*gamma)) - 1);
            df_l = 2 * c_l / (gamma - 1) * (gamma-1)/(2*gamma) * ...
                    (p_guess/p_l)^((gamma-1)/(2*gamma)-1)/p_l;
            
            % 右侧激波函数
            A = 2/((gamma+1)*rho_r);
            B = p_r*(gamma-1)/(gamma+1);
            f_r = (p_guess - p_r) * sqrt(A/(p_guess + B));
            df_r = sqrt(A/(p_guess + B)) * (1 - (p_guess - p_r)/(2*(p_guess + B)));
            
            % 函数总值和导数
            f_total = f_l + f_r + (u_r - u_l);
            df_total = df_l + df_r;
            
            % 更新压力
            p_new = p_guess - f_total / df_total;
            
            % 检查收敛
            if abs(p_new - p_guess) < 1e-6
                p_star = p_new;
                u_star = u_l - f_l;
                return;
            end
            
            p_guess = p_new;
        end
        
        p_star = p_guess;
        u_star = u_l - f_l;
    end
end