function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 3
    T = [T; NaN(3 - size(T, 1), 1)];
end
T(1) = 1/(1+params(1))*1/((1+y(30))*y(28));
T(2) = y(5)^params(4);
T(3) = y(14)^(1-params(4));
end
