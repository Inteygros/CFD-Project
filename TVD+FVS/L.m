function LU = L(U, gamma, dx, FVSfun, TVDLimitersfun)
[~, N] = size(U);
LU = zeros(3, N-2);
[UL, UR] = MUSCL(U, gamma, TVDLimitersfun);
for i = 1:N-2
    [Fplus1,~] = FVSfun(UL(:, i+1), gamma);
    [Fplus2,~] = FVSfun(UL(:, i), gamma);
    [~,Fminus3] = FVSfun(UR(:, i+1), gamma);
    [~,Fminus4] = FVSfun(UR(:, i), gamma);
    LU(:, i) = (Fplus1 + Fminus3 - Fplus2 - Fminus4)/dx;
end
end