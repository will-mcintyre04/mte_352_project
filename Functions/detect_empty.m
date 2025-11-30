function [value, isterminal, direction] = detect_empty(t, h)
    value = h - 1e-3;
    isterminal = 1;
    direction = -1;
end