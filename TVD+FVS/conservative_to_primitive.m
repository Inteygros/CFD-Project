function [rho, u, p, c] = conservative_to_primitive(U, gamma)
    rho = U(1,:);
    u = U(2,:)./U(1,:);
    p = (gamma-1)*(U(3,:)-0.5*U(2,:).^2./U(1,:));
    c = sqrt(gamma*p./U(1,:));
end