% 3½×runge-kutta·½·¨
function U = rk3_step(U, dt, Lfunc)
U0 = U;
L1 = Lfunc(U0, dt);
U(:,2:end-1) = U0(:,2:end-1) - dt * L1;
L2 = Lfunc(U, dt);
U(:,2:end-1) = (3/4) * U0(:,2:end-1) + (1/4) * (U(:,2:end-1) - dt * L2);
L3 = Lfunc(U, dt);
U(:,2:end-1) = (1/3) * U0(:,2:end-1) + (2/3) * (U(:,2:end-1) - dt * L3);
end
