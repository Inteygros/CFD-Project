% Sod���⾫ȷ��
function U = sod_exact_solution(t, N, gamma)
    x = linspace(-0.5, 0.5, N);
    
    % ��ʼ״̬
    rho_l = 1.0;
    u_l = 0.0;
    p_l = 1.0;
    rho_r = 0.125;
    u_r = 0.0;
    p_r = 0.1;
    
    % ���ټ���
    c_l = sqrt(gamma * p_l / rho_l);
    c_r = sqrt(gamma * p_r / rho_r);
    
    % �������������м�״̬
    [p_star, u_star] = solve_riemann(rho_l, u_l, p_l, rho_r, u_r, p_r, gamma);
    
    % �м�״̬�ܶ�
    rho_l_star = rho_l * (p_star/p_l)^(1/gamma); % ���ϡ�貨
    rho_r_star = rho_r * ( (gamma+1)*p_star/p_r + gamma - 1 ) / ...
                ( (gamma+1) + (gamma-1)*p_star/p_r ); % �Ҳ༤��
    
    % �м�״̬����
    c_l_star = c_l * (p_star/p_l)^((gamma-1)/(2*gamma));
    
    % ���ټ���
    s_head = u_l - c_l;                   % ϡ�貨ͷ
    s_tail = u_star - c_l_star;            % ϡ�貨β
    s_shock = u_r + c_r * sqrt(1 + (gamma+1)/(2*gamma)*(p_star/p_r - 1)); % �����ٶ�
    
    % ��ʼ������
    rho = zeros(1, N);
    u = zeros(1, N);
    p = zeros(1, N);
    
    % ����ÿ��λ�õ���㾫ȷ��
    for i = 1:N
        xi = x(i) / t;  % ����������
        
        if xi <= s_head
            % ����1: ��ֹ��
            rho(i) = rho_l;
            u(i) = u_l;
            p(i) = p_l;
            
        elseif xi <= s_tail
            % ����2: ϡ�貨��
            constant = u_l + 2*c_l/(gamma-1);
            u(i) = xi + (gamma-1)/(gamma+1)*(constant - xi);
            c = (gamma-1)/(gamma+1)*(constant - xi);
            rho(i) = rho_l * (c/c_l)^(2/(gamma-1));
            p(i) = p_l * (rho(i)/rho_l)^gamma;
            
        elseif xi <= u_star
            % ����3: �Ӵ�������
            rho(i) = rho_l_star;
            u(i) = u_star;
            p(i) = p_star;
            
        elseif xi <= s_shock
            % ����4: �Ӵ�����Ҳ�
            rho(i) = rho_r_star;
            u(i) = u_star;
            p(i) = p_star;
            
        else
            % ����5: �Ҿ�ֹ��
            rho(i) = rho_r;
            u(i) = u_r;
            p(i) = p_r;
        end
    end
    
    % ת��Ϊ�غ����
    E = p/(gamma-1) + 0.5*rho.*u.^2;
    U = [rho; rho.*u; E];

    % ���������
    function [p_star, u_star] = solve_riemann(rho_l, u_l, p_l, rho_r, u_r, p_r, gamma)
        % ţ�ٵ�����
        max_iter = 20;
        p_guess = 0.3;
        
        for iter = 1:max_iter
            % ���ϡ�貨����
            c_l = sqrt(gamma * p_l / rho_l);
            f_l = 2 * c_l / (gamma - 1) * ((p_guess/p_l)^((gamma-1)/(2*gamma)) - 1);
            df_l = 2 * c_l / (gamma - 1) * (gamma-1)/(2*gamma) * ...
                    (p_guess/p_l)^((gamma-1)/(2*gamma)-1)/p_l;
            
            % �Ҳ༤������
            A = 2/((gamma+1)*rho_r);
            B = p_r*(gamma-1)/(gamma+1);
            f_r = (p_guess - p_r) * sqrt(A/(p_guess + B));
            df_r = sqrt(A/(p_guess + B)) * (1 - (p_guess - p_r)/(2*(p_guess + B)));
            
            % ������ֵ�͵���
            f_total = f_l + f_r + (u_r - u_l);
            df_total = df_l + df_r;
            
            % ����ѹ��
            p_new = p_guess - f_total / df_total;
            
            % �������
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