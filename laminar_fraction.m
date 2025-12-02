function laminar_fraction()
    
    addpath("Functions\");

    % Constants
    g = 9.81;
    rho = 998;
    mu = 0.001003;
    
    d_tube = 0.00525;
    D = 0.110;
    
    epsilon = 1e-6;
    K_entrance = 0.78;
    acceleration_term = 1;
    
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

    laminar_fraction_results = zeros(size(L_values));
    
    h0 = 0.14;

    t_span = [0 300];
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'Events', @detect_empty, 'MaxStep', 0.01);
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
        
        % 1. Find index where it switches from Turbulent to Transition (Re < 4000)
        idx_trans = find(Re < 4000, 1);
        
        % 2. Find index where it switches from Transition to Laminar (Re < 2300)
        idx_laminar = find(Re < 2300, 1);
        
        % Calculate Fraction
        if isempty(idx_laminar)
            laminar_fraction = 0;
        else
            time_laminar = t(end) - t(idx_laminar);
            laminar_fraction = time_laminar / t(end);
            laminar_fraction_results(i) = laminar_fraction;
        end
        
        % Plot the main curve
        plot(t, h, 'LineWidth', 2, 'Color', colors(i,:), ...
             'DisplayName', sprintf('L = %.2f m', p.L));
         
        % --- PLOT MARKERS ---
        % Square for start of Transition (leaving Turbulent)
        if ~isempty(idx_trans)
             plot(t(idx_trans), h(idx_trans), 's', 'MarkerSize', 6, ...
                 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
                 'HandleVisibility', 'off'); 
        end
        
        % Red Circle for start of Laminar (leaving Transition)
        if ~isempty(idx_laminar)
             plot(t(idx_laminar), h(idx_laminar), 'o', 'MarkerSize', 6, ...
                 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
                 'HandleVisibility', 'off'); 
        end
        
        fprintf('     %.2f       |     %.2f     |      %.2f\n', p.L, t(end), laminar_fraction);
    end

    grid on;
    xlabel('Time (seconds)');
    ylabel('Water Height (m)');
    legend('show');
    hold off;
end