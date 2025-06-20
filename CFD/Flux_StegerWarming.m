% Steger-Warming矢通量分裂
function [Fplus, Fminus] = Flux_StegerWarming(U, gamma)
Fplus = zeros(3,1); Fminus = zeros(3,1);
[rho, u, ~, c] = conservative_to_primitive(U, gamma);
lambda1 = u; lambda2 = u+c; lambda3 = u-c;

% 分解成正、负特征值
lambda1p = (lambda1+sqrt(lambda1^2+1e-16))/2; lambda1m = (lambda1-sqrt(lambda1^2+1e-16))/2;
lambda2p = (lambda2+sqrt(lambda2^2+1e-16))/2; lambda2m = (lambda2-sqrt(lambda2^2+1e-16))/2;
lambda3p = (lambda3+sqrt(lambda3^2+1e-16))/2; lambda3m = (lambda3-sqrt(lambda3^2+1e-16))/2;

% 计算正、负通量
Fplus(1) = 2*(gamma-1)*lambda1p+lambda2p+lambda3p;
Fminus(1) = 2*(gamma-1)*lambda1m+lambda2m+lambda3m;
Fplus(2) = 2*(gamma-1)*lambda1p*u+lambda2p*(u+c)+lambda3p*(u-c);
Fminus(2) = 2*(gamma-1)*lambda1m*u+lambda2m*(u+c)+lambda3m*(u-c);
Fplus(3) = (gamma-1)*lambda1p*u^2+(3-gamma)*(lambda2p+lambda3p)*c^2/2/(gamma-1)+lambda2p*(u+c)^2/2+lambda3p*(u-c)^2/2;
Fminus(3) = (gamma-1)*lambda1m*u^2+(3-gamma)*(lambda2m+lambda3m)*c^2/2/(gamma-1)+lambda2m*(u+c)^2/2+lambda3m*(u-c)^2/2;
Fplus = Fplus*rho/2/gamma;
Fminus = Fminus*rho/2/gamma;
end