function [wL, wR] = WENO_5(w0, w1, w2, w3, w4, w5)
wL=zeros(3,1);
wR=zeros(3,1);
for i=1:3
    q0=w0(i)/3-7*w1(i)/6+11*w2(i)/6;
    q1=-w1(i)/6+5*w2(i)/6+w3(i)/3;
    q2=w2(i)/3+5*w3(i)/6-w4(i)/6;
    a0=0.1; a1=0.6; a2=0.3;
    b0=13/12*(w0(i)-2*w1(i)+w2(i))^2+1/4*(w0(i)-4*w1(i)+3*w2(i))^2;
    b1=13/12*(w1(i)-2*w2(i)+w3(i))^2+1/4*(w1(i)-w3(i))^2;
    b2=13/12*(w2(i)-2*w3(i)+w4(i))^2+1/4*(3*w2(i)-4*w3(i)+w4(i))^2;
    omega0=a0/(1e-6+b0)^2;
    omega1=a1/(1e-6+b1)^2;
    omega2=a2/(1e-6+b2)^2;
    omega=omega0+omega1+omega2;
    omega0=omega0/omega;
    omega1=omega1/omega;
    omega2=omega2/omega;
    wL(i)=omega0*q0+omega1*q1+omega2*q2;
    
    q0=w5(i)/3-7*w4(i)/6+11*w3(i)/6;
    q1=-w4(i)/6+5*w3(i)/6+w2(i)/3;
    q2=w3(i)/3+5*w2(i)/6-w1(i)/6;
    a0=0.1; a1=0.6; a2=0.3;
    b0=13/12*(w5(i)-2*w4(i)+w3(i))^2+1/4*(w5(i)-4*w4(i)+3*w3(i))^2;
    b1=13/12*(w4(i)-2*w3(i)+w2(i))^2+1/4*(w4(i)-w2(i))^2;
    b2=13/12*(w3(i)-2*w2(i)+w1(i))^2+1/4*(3*w3(i)-4*w2(i)+w1(i))^2;
    omega0=a0/(1e-6+b0)^2;
    omega1=a1/(1e-6+b1)^2;
    omega2=a2/(1e-6+b2)^2;
    omega=omega0+omega1+omega2;
    omega0=omega0/omega;
    omega1=omega1/omega;
    omega2=omega2/omega;
    wR(i)=omega0*q0+omega1*q1+omega2*q2;
end
end