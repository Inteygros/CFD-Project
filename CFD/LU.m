function LU = LU(U, gamma, dx, FVSfun, FDSfun, TVDLimitersfun, GVCscheme, WENOscheme, UsingCharacteristicReconstruction, Option1, Option2)
[~, N] = size(U);
LU = zeros(3, N-2);

if (Option2==1)
    %使用FVS方法
    if UsingCharacteristicReconstruction %使用特征重构
        % Roe 平均状态
        U_ = Roe_Method(U(:,1:N-1), U(:,2:N), gamma);
        Rcell = cell(1, N-1);
        WL = zeros(3, N-1);
        WR = zeros(3, N-1);
        
        for i = 2:N-2
            % Roe 平均状态的特征基
            [L, R] = Eigen_LR(U_(:, i), gamma);
            
            % 投影到特征空间（Roe 基底）
            W1 = L * U(:, i-1);
            W2 = L * U(:, i);
            W3 = L * U(:, i+1);
            W4 = L * U(:, i+2);
            if i==2
                W0 = W1;
            else
                W0 = L * U(:, i-2);
            end
            if i==N-2
                W5 = W4;
            else
                W5 = L * U(:, i+3);
            end
            if Option1==1
                [WL(:, i), WR(:, i)] = TVD_Format(W1, W2, W3, W4, TVDLimitersfun);
            elseif Option1==2
                [WL(:, i), WR(:, i)] = GVCscheme(W0, W1, W2, W3, W4, W5);
            elseif Option1==3
                [WL(:, i), WR(:, i)] = WENOscheme(W0, W1, W2, W3, W4, W5);
            end
            % 存储R
            Rcell{i} = R;
        end
        
        % 边界点处理（用原始值）
        [~, Rcell{1}] = Eigen_LR(U_(:, 1), gamma);
        [~, Rcell{N-1}] = Eigen_LR(U_(:, N-1), gamma);
        WL(:,1) = WL(:,2);
        WR(:,1) = WR(:,2);
        WR(:,N-1) = WR(:,N-2);
        WL(:,N-1) = WL(:,N-2);
        
        for i = 1:N-2
            [~,~, Aplus1, Aminus1] = FVSfun(U_(:, i+1), gamma);
            [~,~, Aplus2, Aminus2] = FVSfun(U_(:, i), gamma);
            
            R1 = Rcell{i+1};
            R2 = Rcell{i};
            
            Fplus1  = R1 * Aplus1  * WL(:, i+1);
            Fminus1 = R1 * Aminus1 * WR(:, i+1);
            Fplus2  = R2 * Aplus2  * WL(:, i);
            Fminus2 = R2 * Aminus2 * WR(:, i);
            
            LU(:, i) = (Fplus1 + Fminus1 - Fplus2 - Fminus2)/dx;
        end
        
    else
        [UL, UR] = U_Reconstraction(U, TVDLimitersfun, GVCscheme, WENOscheme, Option1);
        for i = 1:N-2
            [Fplus1, ~, ~, ~] = FVSfun(UL(:, i+1), gamma);
            [Fplus2, ~, ~, ~] = FVSfun(UL(:, i), gamma);
            [~, Fminus1, ~, ~] = FVSfun(UR(:, i+1), gamma);
            [~, Fminus2, ~, ~] = FVSfun(UR(:, i), gamma);
            LU(:, i) = (Fplus1 + Fminus1 - Fplus2 - Fminus2)/dx;
        end
    end
else
    %使用FDS方法
    [UL, UR] = U_Reconstraction(U, TVDLimitersfun, GVCscheme, WENOscheme, Option1);
    for i = 1:N-2
        F1 = FDSfun(UL(:, i+1), UR(:, i+1), gamma);
        F2 = FDSfun(UL(:, i), UR(:, i), gamma);
        LU(:, i) = (F1 - F2)/dx;
    end
end
end
