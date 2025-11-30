function V = get_velocity(h, p)
    % Initial guess for velocity with no friction loss
    denom = p.acceleration_term + p.K;
    if (denom == 0)
        denom = 1e-8;
    end
    V = sqrt((2 * p.g * h) / denom);
    
    tol = 1e-8;
    max_iter = 50;
    err = 1;
    iter = 0;
    
    % Iterate to find the correct velocity for a given height
    while err > tol && iter < max_iter
        V_old = V;
        
        Re = (p.rho * V * p.d) / p.mu;
        
        if Re < 2300
            % Laminar
            f = 64 / max(Re, 0.1); 
        elseif Re > 4000
            % Turbulent (Haaland)
            term1 = (p.eps / p.d) / 3.7;
            term2 = 6.9 / Re;
            f = ( -1.8 * log10(term1^1.11 + term2) )^(-2);
        else
            % Transition Zone (2300 <= Re <= 4000)
            % Linear Interpolation logic
            f_2300 = 64 / 2300;
        
            term1 = (p.eps / p.d) / 3.7;
            term2 = 6.9 / 4000;
            f_4000 = ( -1.8 * log10(term1^1.11 + term2) )^(-2);
            
            ratio = (Re - 2300) / (4000 - 2300);
            f = f_2300 + ratio * (f_4000 - f_2300);
        end
        
        denominator = p.K + p.acceleration_term + f * (p.L / p.d);
        V = sqrt( (2 * p.g * h) / denominator );
        
        err = abs(V - V_old);
        iter = iter + 1;
    end
end