% 5阶WENO格式
function f = WENO_5(f0, f1, f2, f3, f4)
f=zeros(3,1);

for i=1:3
    % 插值
    q0=f0(i)/3-7*f1(i)/6+11*f2(i)/6;
    q1=-f1(i)/6+5*f2(i)/6+f3(i)/3;
    q2=f2(i)/3+5*f3(i)/6-f4(i)/6;
    a0=0.1; a1=0.6; a2=0.3;
    
    % 光滑指示器
    b0=13/12*(f0(i)-2*f1(i)+f2(i))^2+1/4*(f0(i)-4*f1(i)+3*f2(i))^2;
    b1=13/12*(f1(i)-2*f2(i)+f3(i))^2+1/4*(f1(i)-f3(i))^2;
    b2=13/12*(f2(i)-2*f3(i)+f4(i))^2+1/4*(3*f2(i)-4*f3(i)+f4(i))^2;
    
    % 权系数
    omega0=a0/(1e-6+b0)^2;
    omega1=a1/(1e-6+b1)^2;
    omega2=a2/(1e-6+b2)^2;
    omega=omega0+omega1+omega2;
    omega0=omega0/omega;
    omega1=omega1/omega;
    omega2=omega2/omega;
    
    f(i)=omega0*q0+omega1*q1+omega2*q2;
end
end