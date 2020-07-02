function [pa, pb] = backward_implicit(zn, pb, pbn, dt, kappa, lambda1, MM, dF)
    pa = zn * (MM + dt*kappa*dF) + lambda1*dt*pbn;
    pb = pb + pbn;
end

