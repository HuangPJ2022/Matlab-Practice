function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 8
    T = [T; NaN(8 - size(T, 1), 1)];
end
T(1) = 1+1/params(7);
T(2) = y(18)-y(14)^T(1)/T(1);
T(3) = y(28)-y(24)^T(1)/T(1);
T(4) = (y(19)*y(5))^params(4);
T(5) = y(14)^(1-params(4));
T(6) = y(1)^params(5)*exp(x(1));
T(7) = (y(19)*y(5))^(params(4)-1);
T(8) = y(14)^(-params(4));
end
