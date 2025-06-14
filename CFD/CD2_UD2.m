function f = CD2_UD2(f1, f2, f3, sig)
f=zeros(3,1);

x=f2-f1; y=f3-f2;
f=f2+sig/4*(sign(x)+sign(y)).*min(abs(x),abs(y));

end