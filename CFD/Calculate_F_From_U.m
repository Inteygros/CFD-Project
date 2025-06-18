% 计算U对应的F
function F = Calculate_F_From_U(U,gamma)
F=zeros(3,1);
F(1)=U(2);
F(2)=U(2)^2/U(1)+(gamma-1)*(U(3)-0.5*U(2)^2/U(1));
F(3)=U(2)/U(1)*(U(3)+(gamma-1)*(U(3)-0.5*U(2)^2/U(1)));
end