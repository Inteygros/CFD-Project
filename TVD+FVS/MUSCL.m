% 使用了特征重构方法
function [UL,UR] = MUSCL(U,gamma,TVDLimitersfun)
[~,N]=size(U);
UL = zeros(3,N-1);
UR = zeros(3,N-1);
U_ = RoeMethod(U(:,1:N-1),U(:,2:N),gamma); %U(:,i)第i+1/2个界面
wL = zeros(3,1); wR = zeros(3,1);
for i = 2:N-2
    [L1, ~] = EigenLR(U_(:,i-1),gamma); [L2, R2] = EigenLR(U_(:,i),gamma); [L3, ~] = EigenLR(U_(:,i+1),gamma);
    w1 = L1*U(:,i-1); w2 = L2*U(:,i); w3 = L2*U(:,i+1); w4 = L3*U(:,i+2);
    
    r1 = (w2(1)-w1(1))/(w3(1)-w2(1)+1e-6); r2 = (w2(2)-w1(2))/(w3(2)-w2(2)+1e-6); r3 = (w2(3)-w1(3))/(w3(3)-w2(3)+1e-6);
    wL(1) = w2(1)+0.5*TVDLimitersfun(r1)*(w3(1)-w2(1));
    wL(2) = w2(2)+0.5*TVDLimitersfun(r2)*(w3(2)-w2(2));
    wL(3) = w2(3)+0.5*TVDLimitersfun(r3)*(w3(3)-w2(3));
    
    r1 = (w3(1)-w2(1))/(w4(1)-w3(1)+1e-6); r2 = (w3(2)-w2(2))/(w4(2)-w3(2)+1e-6); r3 = (w3(3)-w2(3))/(w4(3)-w3(3)+1e-6);
    wR(1) = w3(1)-0.5*TVDLimitersfun(r1)*(w4(1)-w3(1));
    wR(2) = w3(2)-0.5*TVDLimitersfun(r2)*(w4(2)-w3(2));
    wR(3) = w3(3)-0.5*TVDLimitersfun(r3)*(w4(3)-w3(3));
    
    UL(:,i) = R2*wL; UR(:,i) = R2*wR;
end
%处理边界
UL(:,1) = U_(:,1);
UR(:,N-1) = U_(:,N-1);

[L1, ~] = EigenLR(U_(:,N-2),gamma); [L2, R2] = EigenLR(U_(:,N-1),gamma);
w1 = L1*U(:,N-2); w2 = L2*U(:,N-1); w3 = L2*U(:,N);
r1 = (w2(1)-w1(1))/(w3(1)-w2(1)+1e-6); r2 = (w2(2)-w1(2))/(w3(2)-w2(2)+1e-6); r3 = (w2(3)-w1(3))/(w3(3)-w2(3)+1e-6);
wL(1) = w2(1)+0.5*TVDLimitersfun(r1)*(w3(1)-w2(1));
wL(2) = w2(2)+0.5*TVDLimitersfun(r2)*(w3(2)-w2(2));
wL(3) = w2(3)+0.5*TVDLimitersfun(r3)*(w3(3)-w2(3));
UL(:,N-1) = R2*wL;

[L2, R2] = EigenLR(U_(:,1),gamma); [L3, ~] = EigenLR(U_(:,2),gamma);
w2 = L2*U(:,1); w3 = L2*U(:,2); w4 = L3*U(:,3);
r1 = (w3(1)-w2(1))/(w4(1)-w3(1)+1e-6); r2 = (w3(2)-w2(2))/(w4(2)-w3(2)+1e-6); r3 = (w3(3)-w2(3))/(w4(3)-w3(3)+1e-6);
wR(1) = w3(1)-0.5*TVDLimitersfun(r1)*(w4(1)-w3(1));
wR(2) = w3(2)-0.5*TVDLimitersfun(r2)*(w4(2)-w3(2));
wR(3) = w3(3)-0.5*TVDLimitersfun(r3)*(w4(3)-w3(3));
UR(:,1) = R2*wR;

end