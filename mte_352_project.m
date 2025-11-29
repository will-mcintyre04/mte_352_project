function simulate_drainage()
    g = 9.81;
    rho = 1000;
    mu = 0.001;
    
    d_tube = 0.005;
    D_bottle = 0.10;
    
    epsilon = 1e-6;
    K_entrance = 0.78;
    
    p.g = g; p.rho = rho; p.mu = mu;
    p.d = d_tube;
    p.D_bottle = D_bottle;
    p.eps = epsilon;
    p.K_total = 1.78;

    L_values = 0.25:-0.04:0.01;
    
    h0 = 0.14;
    t_span = [0 100];

    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'Events', @detectEmpty);

    figure(1); clf; hold on;
    colors = lines(length(L_values));
    
    fprintf('Tube Length (m) | Total Time (s)\n');
    fprintf('-------------------------------\n');

    for i = 1:length(L_values)
        p.L = L_values(i);
        
        [t, h] = ode45(@(t,h) drain_physics(t, h, p), t_span, h0, options);
        
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

function dh_dt = drain_physics(t, h, p)
    h_calc = max(0, h); 

    V = sqrt((2 * p.g * h_calc) / p.K_total);
    
    tol = 1e-8;
    max_iter = 50;
    err = 1;
    iter = 0;
    
    while err > tol && iter < max_iter
        V_old = V;
        
        if V < 1e-9
            V = 0;
            f = 0.03;
        else
            Re = (p.rho * V * p.d) / p.mu;
            
            if Re < 2300
                f = 64 / max(Re, 0.1); 
            else
                term1 = (p.eps / p.d) / 3.7;
                term2 = 6.9 / Re;
                f = ( -1.8 * log10(term1^1.11 + term2) )^(-2);
            end
        end
        
        denominator = p.K_total + f * (p.L / p.d);
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