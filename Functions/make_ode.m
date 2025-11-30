function dh_dt = make_ode(h, p)
    h = max(0, h); 
    V = get_velocity(h, p);
    dh_dt = -(p.d / p.D)^2 * V;
end