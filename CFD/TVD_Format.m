function f = TVD_Format(f1, f2, f3, TVDLimitersfun, sig)
f = zeros(3,1);
small = 1e-6;
r = (f2 - f1) ./ (f3 - f2 + small);
phi = TVDLimitersfun(r);
f = f2 + (sig / 2) * phi .* (f3 - f2);
end