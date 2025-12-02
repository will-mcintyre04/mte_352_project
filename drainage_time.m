function drainage_time()
    
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
    
    h0 = 0.14;
    t_span = [0 300];
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'Events', @detect_empty);

    drainage_times = zeros(size(L_values));
    
    for i = 1:length(L_values)
        p.L = L_values(i);
        
        [t, h] = ode45(@(t, h) make_ode(h, p), t_span, h0, options);
        
        drainage_times(i) = t(end);
    end

    figure(2); clf;
    set(gcf, 'Color', "white", "Name", "Drainage Time");
    plot(L_values.*100, drainage_times, "Marker",".","MarkerSize",20);
    grid on;
    ylim([0 max(drainage_times)+10]);
    xlabel('Straw Length [cm]');
    ylabel('Drainage Time [seconds]');
    hold off;
end