function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(5, 1);
end
[T_order, T] = HPJedit_baseline_RBC.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(10, 1);
    residual(1) = (T(1)) - (T(1)*params(1)*(1+y(2)*y(9)-y(10)));
    residual(2) = (y(4)) - ((y(3)/y(8))^params(7));
    residual(3) = (y(7)) - (y(1)*T(2)*T(3));
    residual(4) = (y(7)) - (y(8)+y(6));
    residual(5) = (y(5)) - (y(6)+y(5)*(1-y(10)));
    residual(6) = (y(1)) - (y(1)^params(5)*exp(x(1)));
    residual(7) = (y(2)) - (T(3)*y(1)*params(4)*T(4));
    residual(8) = (y(2)) - (params(2)*params(6)*y(9)^(params(6)-1));
    residual(9) = (y(3)) - (T(2)*y(1)*(1-params(4))*T(5));
    residual(10) = (y(10)) - (params(2)*y(9)^params(6));
end
