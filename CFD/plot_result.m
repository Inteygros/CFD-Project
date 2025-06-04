function plot_result(U, h, gamma)
    [rho, u, p, ~] = conservative_to_primitive(U, gamma);
    set(h(1), 'YData', rho);
    set(h(2), 'YData', u);
    set(h(3), 'YData', p);
end