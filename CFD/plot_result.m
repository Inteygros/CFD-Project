% 绘图函数
function plot_result(U, h, gamma, t, N, show_exact_solution)
    % 更新数值解
    [rho, u, p, ~] = conservative_to_primitive(U, gamma);
    set(h{1,2}, 'YData', rho);
    set(h{2,2}, 'YData', u);
    set(h{3,2}, 'YData', p);
    
    % 计算并绘制精确解
    if show_exact_solution
    exact_U = sod_exact_solution(t, N, gamma);
    [exact_rho, exact_u, exact_p, ~] = conservative_to_primitive(exact_U, gamma);
    set(h{1,1}, 'YData', exact_rho);
    set(h{2,1}, 'YData', exact_u);
    set(h{3,1}, 'YData', exact_p);
    end
    
    % 更新标题显示当前时间
    title(sprintf('Sod Shock Tube Solution at t = %.4f', t));
end