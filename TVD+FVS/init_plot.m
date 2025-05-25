function h = init_plot(x, rho, u, p)
    h(1) = plot(x, rho, 'r-', 'LineWidth', 1.5); hold on;
    h(2) = plot(x, u,   'g-', 'LineWidth', 1.5);
    h(3) = plot(x, p,   'b-', 'LineWidth', 1.5);
    xlim([-0.5, 0.5]);
    ylim([0, 1]);
    xlabel('x');
    ylabel('Value');
    legend(h, {'\rho (Density)', 'u (Velocity)', 'p (Pressure)'}, 'Location', 'northeast');
    title('Sod Shock Tube Solution');
    grid on;
end