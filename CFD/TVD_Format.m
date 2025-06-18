% TVD格式
function f = TVD_Format(f1, f2, f3, TVDLimitersfun, c, sig)
f = zeros(3,1);
small = 1e-10;
if sig>0 % 正通量
    r = (f2 - f1) ./ (f3 - f2 + small);
    phi = TVDLimitersfun(r);
    f = f2 + 1 / 2 * phi .* (1-c) .*(f3 - f2);
else % 负通量
    r = (f3 - f2) ./ (f2 - f1 + small);
    phi = TVDLimitersfun(r);
    f = f2 - 1 / 2 * phi .* (1+c) .*(f2 - f1);
end
end