% rows 25 -> 1.3 cm straw length
% colums 14 -> 2 cm water height
t = [
    0	8.79	19.61	31.45	44.4	59.68	80.41 ;
    0	8.6	    18.39	30.8	44.53	58.55	81.63 ;
    0	8.07	18.84	30.32	43.25	57.48	77.74 ;
    0	8.77	18.3	29.35	42.82	57.09	78.89 ;
    0	7.14	17.46	26.38	38.65	52.38	70.99 ;
    0	7.574	15.749	24.691	36.103	48.148	67.467;
    0	6.263	14.648	23.556	34.234	45.945	63.923;
];
z1s = (14:-2:2)./100;
Ls  = (25:-4:1)./100;

colors = lines(length(Ls));

drainage_profiles();

figure(1);
hold on;
for i=1:length(Ls)
    plot(t(i,:), z1s, "Color", colors(i,:), "LineWidth", 1.5, ...
        'DisplayName', sprintf('L = %.2f m', Ls(i)), ...
             'LineStyle','--');
end
grid on;
title("")
xlabel('Time [seconds]');
ylabel('Water Height [m]');
legend('show');
hold off;