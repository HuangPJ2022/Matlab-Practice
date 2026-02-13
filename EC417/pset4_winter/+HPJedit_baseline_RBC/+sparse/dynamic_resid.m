function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(6, 1);
end
[T_order, T] = HPJedit_baseline_RBC.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(10, 1);
    residual(1) = (1/y(18)) - (params(1)*(1+y(22)*y(29)-y(30))*T(1));
    residual(2) = (y(14)) - ((y(13)/y(18))^params(7));
    residual(3) = (y(17)) - (y(11)*T(2)*T(3));
    residual(4) = (y(17)) - (y(18)+y(16));
    residual(5) = (y(15)) - (y(16)+y(5)*(1-y(20)));
    residual(6) = (y(11)) - (T(4));
    residual(7) = (y(12)) - (T(3)*y(11)*params(4)*T(5));
    residual(8) = (y(12)) - (params(2)*params(6)*y(19)^(params(6)-1));
    residual(9) = (y(13)) - (T(2)*y(11)*(1-params(4))*T(6));
    residual(10) = (y(20)) - (params(2)*y(19)^params(6));
end
