function [wL,wR] = CD2_UD2(~, w1, w2, w3, w4, ~)
wL=zeros(3,1);
wR=zeros(3,1);
for i=1:3
    if abs(w1(i)-w2(i))<abs(w2(i)-w3(i))
        wL(i)=0.5*(3*w2(i)-w1(i));
    else
        wL(i)=0.5*(w2(i)+w1(i));
    end
    if abs(w2(i)-w3(i))<abs(w3(i)-w4(i))
        wR(i)=0.5*(3*w3(i)-w2(i));
    else
        wR(i)=0.5*(w3(i)+w2(i));
    end
end
end