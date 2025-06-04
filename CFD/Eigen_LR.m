function [L, R] = Eigen_LR(U, gamma)
[rho, u, p, c] = conservative_to_primitive(U, gamma);
E = U(3);
H = (E + p) / rho;

% 右特征向量矩阵 R（按特征值顺序 u, u+c, u-c）
R = [1, 1, 1;
    u, u + c, u - c;
    0.5*u^2, H + u*c, H - u*c];

c2 = c^2;
% 左特征向量矩阵 R
L = zeros(3,3);

L(1,1) = 1 - (gamma-1)*u^2/(2*c2);
L(1,2) = (gamma-1)*u/c2;
L(1,3) = -(gamma-1)/c2;

L(2,1) = (gamma-1)*u^2/(4*c2) - u/(2*c);
L(2,2) = -(gamma-1)*u/(2*c2) + 1/(2*c);
L(2,3) = (gamma-1)/(2*c2);

L(3,1) = (gamma-1)*u^2/(4*c2) + u/(2*c);
L(3,2) = -(gamma-1)*u/(2*c2) - 1/(2*c);
L(3,3) = (gamma-1)/(2*c2);

end