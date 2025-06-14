function h = init_plot(x, U, gamma)
[rho, u, p, ~] = conservative_to_primitive(U, gamma);
h = cell(3, 2);
h{1,1} = plot(x, rho, 'r--', 'LineWidth', 0.8); hold on;  % 精确解 - 虚线
h{1,2} = plot(x, rho, 'r-', 'LineWidth', 0.8);          % 数值解 - 实线
h{2,1} = plot(x, u, 'g--', 'LineWidth', 0.8);
h{2,2} = plot(x, u, 'g-', 'LineWidth', 0.8);
h{3,1} = plot(x, p, 'b--', 'LineWidth', 0.8);
h{3,2} = plot(x, p, 'b-', 'LineWidth', 0.8);
xlim([-0.5, 0.5]);
ylim([-0.1, 1.1]);
xlabel('x');
ylabel('Value');
title('Sod Shock Tube Solution');
grid on;
legend([h{1,1}, h{1,2}, h{2,1}, h{2,2}, h{3,1}, h{3,2}], ...
           {'Exact \rho', 'Numerical \rho', 'Exact u', 'Numerical u', ...
            'Exact p', 'Numerical p'}, 'Location', 'northeast');

end