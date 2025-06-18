% MCÏÞÖÆÆ÷
function phi = MC_limiter(r)
    phi = max(0,min(min((1+r)/2,2),2*r));
end