function [L, R] = EigenLR(U, gamma)
[rho, u, p, c] = conservative_to_primitive(U, gamma);
E = U(3);
H = (E + p) / rho;

% 右特征向量矩阵 R（按特征值顺序 u, u+c, u-c）
R = [1, 1, 1;
    u, u + c, u - c;
    0.5*u^2, H + u*c, H - u*c];

% 左特征向量矩阵 R
den = 2*H - u^2;
L = zeros(3,3);
L(1,1) = 2*(H - u^2)/den;
L(1,2) = 2*u/den;
L(1,3) = 2/(u^2 - 2*H);

L(2,1) = u*(-H + 0.5*c*u + 0.5*u^2)/(c * den);
L(2,2) = (H - c*u - 0.5*u^2)/(c * den);
L(2,3) = 1 / den;

L(3,1) = u*(H + 0.5*c*u - 0.5*u^2)/(c * den);
L(3,2) = (-H - c*u + 0.5*u^2)/(c * den);
L(3,3) = 1 / den;

end