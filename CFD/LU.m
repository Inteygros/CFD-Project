function LU = LU(U, gamma, dx, dt, FVSfun, FDSfun, TVDLimitersfun, GVCscheme, WENOscheme, UsingCharacteristicReconstruction, Option1, Option2)
[~, N] = size(U);
LU = zeros(3, N-2);
Flux = zeros(3, N-1);
% Roe 平均状态
U_ = Roe_Method(U(:,1:N-1), U(:,2:N), gamma);

if Option2==1
    %使用FVS方法
    Fplus=zeros(3,N);
    Fminus=zeros(3,N);
    for i = 1:N
        [Fplus(:, i), Fminus(:, i)] = FVSfun(U(:, i), gamma);
    end
    FL = zeros(3, N-1);
    FR = zeros(3, N-1);
    
    if UsingCharacteristicReconstruction %使用特征重构

        Rcell = cell(1, N-1);
        
        for i = 2:N-2
            % Roe 平均状态的特征基
            [L, R] = Eigen_LR(U_(:, i), gamma);
            
            % 投影到特征空间
            f1plus = L * Fplus(:, i-1);
            f1minus = L * Fminus(:, i-1);
            f2plus = L * Fplus(:, i);
            f2minus = L * Fminus(:, i);
            f3plus = L * Fplus(:, i+1);
            f3minus = L * Fminus(:, i+1);
            f4plus = L * Fplus(:, i+2);
            f4minus = L * Fminus(:, i+2);
            if i==2
                f0plus = f1plus;
            else
                f0plus = L * Fplus(:,i-2);
            end
            if i==N-2
                f5minus = f4minus;
            else
                f5minus = L * Fminus(:,i+3);
            end
            if Option1==1 % TVD
                [~, u, ~, c] = conservative_to_primitive(U_(:,i), gamma);
                lambda1 = u; lambda2 = u+c; lambda3 = u-c;
                cp = dt/dx*[(lambda1+sqrt(lambda1^2+1e-10))/2;
                    (lambda2+sqrt(lambda2^2+1e-10))/2;
                    (lambda3+sqrt(lambda3^2+1e-10))/2];
                cm = dt/dx*[(lambda1-sqrt(lambda1^2+1e-10))/2;
                    (lambda2-sqrt(lambda2^2+1e-10))/2;
                    (lambda3-sqrt(lambda3^2+1e-10))/2];
                FL(:,i) = TVD_Format(f1plus, f2plus, f3plus, TVDLimitersfun, cp, 1);
                FR(:,i) = TVD_Format(f2minus, f3minus, f4minus, TVDLimitersfun, cm, -1);
            elseif Option1==2 % GVC
                FL(:,i) = GVCscheme(f1plus, f2plus, f3plus, 1);
                FR(:,i) = GVCscheme(f2minus, f3minus, f4minus, -1);
            elseif Option1==3 % WENO
                FL(:,i) = WENOscheme(f0plus, f1plus, f2plus, f3plus, f4plus);
                FR(:,i) = WENOscheme(f5minus, f4minus, f3minus, f2minus, f1minus);
            end
            % 存储R
            Rcell{i} = R;
        end
        
        % 边界点处理（用原始值）
        [~, Rcell{1}] = Eigen_LR(U_(:, 1), gamma);
        [~, Rcell{N-1}] = Eigen_LR(U_(:, N-1), gamma);
        FL(:,1) = FL(:,2);
        FR(:,1) = FR(:,2);
        FR(:,N-1) = FR(:,N-2);
        FL(:,N-1) = FL(:,N-2);
        
        % 计算LU
        for i = 1:N-1
            Flux(:,i) = Rcell{i}*(FL(:,i)+FR(:,i));
        end
        for i = 1:N-2
            LU(:, i) = (Flux(:,i+1) - Flux(:,i))/dx;
        end
        
    else %不使用特征重构
        for i = 2:N-2
            if i==2
                f0plus = Fplus(:,i-1);
            else
                f0plus = Fplus(:,i-2);
            end
            if i==N-2
                f5minus = Fminus(:,i+2);
            else
                f5minus = Fminus(:,i+3);
            end
            if Option1==1 % TVD
                [~, u, ~, c] = conservative_to_primitive(U_(:,i), gamma);
                lambda1 = u; lambda2 = u+c; lambda3 = u-c;
                cp = dt/dx*[(lambda1+sqrt(lambda1^2+1e-10))/2;
                    (lambda2+sqrt(lambda2^2+1e-10))/2;
                    (lambda3+sqrt(lambda3^2+1e-10))/2];
                cm = dt/dx*[(lambda1-sqrt(lambda1^2+1e-10))/2;
                    (lambda2-sqrt(lambda2^2+1e-10))/2;
                    (lambda3-sqrt(lambda3^2+1e-10))/2];
                FL(:,i) = TVD_Format(Fplus(:,i-1), Fplus(:,i), Fplus(:,i+1), TVDLimitersfun, cp, 1);
                FR(:,i) = TVD_Format(Fminus(:,i), Fminus(:,i+1), Fminus(:,i+2), TVDLimitersfun, cm, -1);
            elseif Option1==2 % GVC
                FL(:,i) = GVCscheme(Fplus(:,i-1), Fplus(:,i), Fplus(:,i+1), 1);
                FR(:,i) = GVCscheme(Fminus(:,i), Fminus(:,i+1), Fminus(:,i+2), -1);
            elseif Option1==3 % WENO
                FL(:,i) = WENOscheme(f0plus, Fplus(:,i-1), Fplus(:,i), Fplus(:,i+1), Fplus(:,i+2));
                FR(:,i) = WENOscheme(f5minus, Fminus(:,i+2), Fminus(:,i+1), Fminus(:,i), Fminus(:,i-1));
            end
        end
        %处理边界
        FL(:,1)=FL(:,2);
        FR(:,1)=FR(:,2);
        FL(:,N-1)=FL(:,N-2);
        FR(:,N-1)=FR(:,N-2);
        Flux=FL+FR;
    end
    for i = 1:N-2
        LU(:, i) = (Flux(:,i+1) - Flux(:,i))/dx;
    end
end

if Option2==2
    %使用FDS方法
    UL = zeros(3, N-1);
    UR = zeros(3, N-1);
    for i = 2:N-2
        if i==2
            U0 = U(:,i-1);
        else
            U0 = U(:,i-2);
        end
        if i==N-2
            U5 = U(:,i+2);
        else
            U5 = U(:,i+3);
        end
        if Option1==1 % TVD
            [~, u, ~, c] = conservative_to_primitive(U_(:,i), gamma);
                lambda1 = u; lambda2 = u+c; lambda3 = u-c;
                cp = dt/dx*[(lambda1+sqrt(lambda1^2+1e-6))/2;
                    (lambda2+sqrt(lambda2^2+1e-6))/2;
                    (lambda3+sqrt(lambda3^2+1e-6))/2];
                cm = dt/dx*[(lambda1-sqrt(lambda1^2+1e-6))/2;
                    (lambda2-sqrt(lambda2^2+1e-6))/2;
                    (lambda3-sqrt(lambda3^2+1e-6))/2];
            UL(:,i) = TVD_Format(U(:,i-1), U(:,i), U(:,i+1), TVDLimitersfun,cp, 1);
            UR(:,i) = TVD_Format(U(:,i), U(:,i+1), U(:,i+2), TVDLimitersfun,cm, -1);
        elseif Option1==2 % GVC
            UL(:,i) = GVCscheme(U(:,i-1), U(:,i), U(:,i+1), 1);
            UR(:,i) = GVCscheme(U(:,i), U(:,i+1), U(:,i+2), -1);
        elseif Option1==3 % WENO
            UL(:,i) = WENOscheme(U0, U(:,i-1), U(:,i), U(:,i+1), U(:,i+2));
            UR(:,i) = WENOscheme(U5, U(:,i+2), U(:,i+1), U(:,i), U(:,i-1));
        end
    end
    UL(:,1)=UL(:,2);
    UR(:,1)=UR(:,2);
    UL(:,N-1)=UL(:,N-2);
    UR(:,N-1)=UR(:,N-2);
    for i = 1:N-1
        Flux(:, i) = FDSfun(UL(:, i), UR(:, i), gamma);
    end
    
    for i = 1:N-2
        LU(:, i) = (Flux(:,i+1) - Flux(:,i))/dx;
    end
end
end
