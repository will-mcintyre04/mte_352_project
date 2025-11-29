function mte_352_project()

    % Constants
    g = 9.81;
    rho = 998;
    mu = 0.001;
    
    d_tube = 0.00525;
    D_bottle = 0.110;
    
    epsilon = 1e-6;
    K_entrance = 0.78;
    acceleration_term = 1;
    
    % Structure for drain function
    p.g = g; p.rho = rho; p.mu = mu;
    p.d = d_tube;
    p.D_bottle = D_bottle;
    p.eps = epsilon;
    p.K = K_entrance;
    p.acceleration_term = acceleration_term;

    L_values = 0.25:-0.04:0.01;
    
    h0 = 0.14;
    t_span = [0 120];

    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'Events', @detectEmpty);

    figure(1); clf; hold on;
    colors = lines(length(L_values));
    
    fprintf('Tube Length (m) | Total Time (s)\n');
    fprintf('-------------------------------\n');

    for i = 1:length(L_values)
        p.L = L_values(i);
        
        [t, h] = ode45(@(t, h) make_ode(h, p), t_span, h0, options);
        
        plot(t, h, 'LineWidth', 2, 'Color', colors(i,:), ...
             'DisplayName', sprintf('L = %.2f m', p.L));
        
        fprintf('     %.2f       |     %.2f\n', p.L, t(end));
    end

    grid on;
    xlabel('Time (seconds)');
    ylabel('Water Height (m)');
    title('Bottle Drainage Profiles for Various Tube Lengths');
    legend('show');
    hold off;
end

function dh_dt = make_ode(h, p)
    h_calc = max(0, h); 

    % Initial guess for velocity with no friction loss
    V = sqrt((2 * p.g * h_calc) / (p.acceleration_term + p.K));
    
    tol = 1e-8;
    max_iter = 50;
    err = 1;
    iter = 0;
    
    % Iterate to find the correct velocity for a given height
    while err > tol && iter < max_iter
        V_old = V;
        
        Re = (p.rho * V * p.d) / p.mu;
        
        if Re < 2300
            f = 64 / max(Re, 0.1); 
        else
            term1 = (p.eps / p.d) / 3.7;
            term2 = 6.9 / Re;
            f = ( -1.8 * log10(term1^1.11 + term2) )^(-2);
        end
        
        denominator = (p.K + p.acceleration_term) + f * (p.L / p.d);
        V = sqrt( (2 * p.g * h_calc) / denominator );
        
        err = abs(V - V_old);
        iter = iter + 1;
    end
    
    A_tube = pi * (p.d / 2)^2;
    A_bottle = pi * (p.D_bottle / 2)^2;
    
    dh_dt = -(A_tube / A_bottle) * V;
end

function [value, isterminal, direction] = detectEmpty(t, h)
    value = h;
    isterminal = 1;
    direction = -1;
end