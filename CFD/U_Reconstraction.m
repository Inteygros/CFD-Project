function [UL,UR] = U_Reconstraction(U, TVDLimitersfun, GVCscheme, WENOscheme, Option1)
[~,N]=size(U);
UL = zeros(3,N-1);
UR = zeros(3,N-1);

for i = 2:N-2
    w1 = U(:,i-1); w2 = U(:,i); w3 = U(:,i+1); w4 = U(:,i+2);
    if i==2
        w0 = w1;
    else
        w0 = U(:, i-2);
    end
    if i==N-2
        w5 = w4;
    else
        w5 = U(:, i+3);
    end
    if Option1==1 % TVD
        [UL(:,i),UR(:,i)] = TVD_Format(w1, w2, w3, w4, TVDLimitersfun);
    elseif Option1==2 % GVC
        [UL(:,i),UR(:,i)] = GVCscheme(w0, w1, w2, w3, w4, w5);
    elseif Option1==3 % WENO
        [UL(:,i),UR(:,i)] = WENOscheme(w0, w1, w2, w3, w4, w5);
    end
end
%¥¶¿Ì±ﬂΩÁ
UL(:,1)=UL(:,2);
UR(:,1)=UR(:,2);
UL(:,N-1)=UL(:,N-2);
UR(:,N-1)=UR(:,N-2);
end