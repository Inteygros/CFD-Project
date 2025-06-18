% 绘图的初始化
function h = init_plot(x, U, gamma, show_exact_solution)
[rho, u, p, ~] = conservative_to_primitive(U, gamma);
h = cell(3, 2);
h{1,2} = plot(x, rho, 'r-', 'LineWidth', 0.8); hold on;       % 数值解 - 实线
h{2,2} = plot(x, u, 'g-', 'LineWidth', 0.8);
h{3,2} = plot(x, p, 'b-', 'LineWidth', 0.8);

if show_exact_solution
    h{1,1} = plot(x, rho, 'r--', 'LineWidth', 0.5);   % 精确解 - 虚线
    h{2,1} = plot(x, u, 'g--', 'LineWidth', 0.5);
    h{3,1} = plot(x, p, 'b--', 'LineWidth', 0.5);
end

xlim([-0.5, 0.5]);
ylim([-0.1, 1.1]);
xlabel('x');
ylabel('Value');
title('Sod Shock Tube Solution');
grid on;
if show_exact_solution
    legend([h{1,1}, h{1,2}, h{2,1}, h{2,2}, h{3,1}, h{3,2}], ...
        {'Exact \rho', 'Numerical \rho', 'Exact u', 'Numerical u', 'Exact p', 'Numerical p'}, ...
        'Location', 'northeast');
else
    legend([h{1,2}, h{2,2}, h{3,2}], ...
        {'Numerical \rho', 'Numerical u', 'Numerical p'}, ...
        'Location', 'northeast');
end
end