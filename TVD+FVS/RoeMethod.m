function U_ = RoeMethod(UL, UR, gamma)
[~,N]=size(UL);
U_ = zeros(3,N);
for i = 1:N
    [rhoL, uL, pL, ~] = conservative_to_primitive(UL(:,i), gamma);
    [rhoR, uR, pR, ~] = conservative_to_primitive(UR(:,i), gamma);
    cpower2L = gamma*pL/rhoL; cpower2R = gamma*pR/rhoR; 
    HL = cpower2L/(gamma-1)+0.5*uL^2; HR = cpower2R/(gamma-1)+0.5*uR^2;
    rho_ = sqrt(rhoL*rhoR);
    u_ = (sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR));
    H_ = (sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR));
    U_(1,i) = rho_;
    U_(2,i) = rho_*u_;
    U_(3,i) = rho_*(0.5*u_^2*(1-1/gamma)+H_/gamma);
end
end