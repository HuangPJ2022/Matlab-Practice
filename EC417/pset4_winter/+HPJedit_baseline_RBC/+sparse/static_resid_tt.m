function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 5
    T = [T; NaN(5 - size(T, 1), 1)];
end
T(1) = 1/y(8);
T(2) = (y(9)*y(5))^params(4);
T(3) = y(4)^(1-params(4));
T(4) = (y(9)*y(5))^(params(4)-1);
T(5) = y(4)^(-params(4));
end
