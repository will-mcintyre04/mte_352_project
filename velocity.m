addpath("Functions\")

% Constants
g = 9.81;
rho = 998;
mu = 0.001002;

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
t_span = [0 120];
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'Events', @detect_empty);
figure(1); clf; hold on;
set(gcf, 'Color', "white", "Name", "Drainage Profiles");
colors = lines(length(L_values));

fprintf('Tube Length (m) | Total Time (s)\n');
fprintf('--------------------------------\n');
for i = 1:length(L_values)
    p.L = L_values(i);
    
    [t, h] = ode45(@(t, h) make_ode(h, p), t_span, h0, options);
    
    V = get_velocity(h, p);
    % Plot the main curve
    plot(t, h, 'LineWidth', 2, 'Color', colors(i,:), ...
         'DisplayName', sprintf('L = %.2f m', p.L));
    
    fprintf('     %.2f       |     %.2f\n', p.L, t(end));
end
grid on;
xlabel('Time (seconds)');
ylabel('Water Height (m)');
title('Drainage Profiles');
legend('show');
hold off;

