function mte_352_project()

    % Constants
    g = 9.81;
    rho = 998;
    mu = 0.001002;
    
    d_tube = 0.00525;
    D = 0.110;
    
    epsilon = 1e-6;
    K_entrance = 0;
    acceleration_term = 0;
    
    % Structure for drain function
    p.g = g;
    p.rho = rho;
    p.mu = mu;
    p.d = d_tube;
    p.D = D;
    p.eps = epsilon;
    p.K = K_entrance;
    p.acceleration_term = acceleration_term;

    L_values = 0.25:-0.04:0.01;
    
    h0 = 0.14;
    t_span = [0 120];

    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'Events', @detectEmpty);

    figure(1); clf; hold on;
    colors = lines(length(L_values));
    
    fprintf('Tube Length (m) | Total Time (s) | Laminar Fraction\n');
    fprintf('---------------------------------------------------\n');

    for i = 1:length(L_values)
        p.L = L_values(i);
        
        [t, h] = ode45(@(t, h) make_ode(h, p), t_span, h0, options);

        % Calculate all reynolds numbers
        Re = zeros(size(h));
        for k = 1:length(h)
            Re(k) = calculate_Re(h(k), p);
        end

        % find index where it switches from turbulent to laminar(Re <2300)
        idx_laminar = find(Re < 2300, 1);

        time_laminar = t(end) - t(idx_laminar);

        laminar_fraction = time_laminar / t(end);

        plot(t, h, 'LineWidth', 2, 'Color', colors(i,:), ...
             'DisplayName', sprintf('L = %.2f m', p.L));
        
        fprintf('     %.2f       |     %.2f     |      %.2f\n', p.L, t(end), laminar_fraction);
    end

    grid on;
    xlabel('Time (seconds)');
    ylabel('Water Height (m)');
    title('Bottle Drainage Profiles for Various Tube Lengths');
    legend('show');
    hold off;
end

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
            f = 64 / max(Re, 0.1); 
        else
            term1 = (p.eps / p.d) / 3.7;
            term2 = 6.9 / Re;
            f = ( -1.8 * log10(term1^1.11 + term2) )^(-2);
        end
        
        denominator = p.K + p.acceleration_term + f * (p.L / p.d);
        V = sqrt( (2 * p.g * h) / denominator );
        
        err = abs(V - V_old);
        iter = iter + 1;
    end
end

function Re = calculate_Re(h, p)
    h = max(0, h);
    V = get_velocity(h, p);
    Re = p.rho * V * p.d / p.mu;
end

function dh_dt = make_ode(h, p)
    h = max(0, h); 
    V = get_velocity(h, p);
    dh_dt = -(p.d / p.D)^2 * V;
end

function [value, isterminal, direction] = detectEmpty(t, h)
    value = h - 1e-3;
    isterminal = 1;
    direction = -1;
end