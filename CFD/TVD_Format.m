function [wL,wR] = TVD_Format(w1, w2, w3, w4, TVDLimitersfun)
wL=zeros(3,1);
wR=zeros(3,1);
r1 = (w2(1)-w1(1))/(w3(1)-w2(1)+1e-6); r2 = (w2(2)-w1(2))/(w3(2)-w2(2)+1e-6); r3 = (w2(3)-w1(3))/(w3(3)-w2(3)+1e-6);
wL(1) = w2(1)+0.5*TVDLimitersfun(r1)*(w3(1)-w2(1));
wL(2) = w2(2)+0.5*TVDLimitersfun(r2)*(w3(2)-w2(2));
wL(3) = w2(3)+0.5*TVDLimitersfun(r3)*(w3(3)-w2(3));

r1 = (w3(1)-w2(1))/(w4(1)-w3(1)+1e-6); r2 = (w3(2)-w2(2))/(w4(2)-w3(2)+1e-6); r3 = (w3(3)-w2(3))/(w4(3)-w3(3)+1e-6);
wR(1) = w3(1)-0.5*TVDLimitersfun(r1)*(w4(1)-w3(1));
wR(2) = w3(2)-0.5*TVDLimitersfun(r2)*(w4(2)-w3(2));
wR(3) = w3(3)-0.5*TVDLimitersfun(r3)*(w4(3)-w3(3));
end