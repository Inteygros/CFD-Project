function [wL,wR] = UC3_SC4(w0, w1, w2, w3, w4, ~)
% 群速度控制格式（UC3 + SC4）
wL = zeros(3,1);
wR = zeros(3,1);
alpha = 0.7;

for k = 1:3
    % UC3 插值
    wL_UC3 = ( 2*w0(k) - 7*w1(k) + 11*w2(k) ) / 6;
    wR_UC3 = ( 11*w3(k) - 7*w2(k) + 2*w1(k) ) / 6;

    % SC4 插值
    wL_SC4 = (-w0(k) + 7*w1(k) + 7*w2(k) - w3(k)) / 12;
    wR_SC4 = (-w1(k) + 7*w2(k) + 7*w3(k) - w4(k)) / 12;

    % 光滑度指示
    betaL = abs(w0(k)-w1(k)) < abs(w1(k)-w2(k));
    betaR = abs(w1(k)-w2(k)) < abs(w2(k)-w3(k));

    % 插值混合
    wL(k) = betaL * (alpha * wL_UC3 + (1-alpha) * wL_SC4) + (1 - betaL) * wL_SC4;
    wR(k) = betaR * (alpha * wR_UC3 + (1-alpha) * wR_SC4) + (1 - betaR) * wR_SC4;
end
end