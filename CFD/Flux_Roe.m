function flux = Flux_Roe(UL, UR, gamma)
FL = Calculate_F_From_U(UL, gamma);
FR = Calculate_F_From_U(UR, gamma);
U = Roe_Method(UL, UR, gamma);
[L, R] = Eigen_LR(U, gamma);
[~, u, ~, c] = conservative_to_primitive(U, gamma);
abs1=abs(u); abs2=abs(u+c); abs3=abs(u-c);
% 使用了Harten熵修正技术
delta=0.1*c;
if abs1<delta
    abs1=(abs1^2+delta^2)/2/delta;
end
if abs2<delta
    abs2=(abs2^2+delta^2)/2/delta;
end
if abs3<delta
    abs3=(abs3^2+delta^2)/2/delta;
end
absA = R*diag([abs1, abs2, abs3])*L;
flux = 0.5*(FL+FR-absA*(UR-UL));
end