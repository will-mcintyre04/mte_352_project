function Re = calculate_Re(h, p)
    h = max(0, h);
    V = get_velocity(h, p);
    Re = p.rho * V * p.d / p.mu;
end