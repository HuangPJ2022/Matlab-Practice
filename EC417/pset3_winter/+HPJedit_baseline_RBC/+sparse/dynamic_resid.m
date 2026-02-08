function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(3, 1);
end
[T_order, T] = HPJedit_baseline_RBC.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(10, 1);
    residual(1) = (1/((1+y(20))*y(18))) - (T(1)*(1+y(22)-params(2)));
    residual(2) = (params(6)/(1-y(14))) - (y(13)/((1+y(20))*y(18)));
    residual(3) = (y(17)) - (y(11)*T(2)*T(3));
    residual(4) = (y(13)) - (y(17)*(1-params(4))/y(14));
    residual(5) = (y(12)) - (y(17)*params(4)/y(5));
    residual(6) = (y(15)) - (y(5)*(1-params(2))+y(16));
    residual(7) = (y(17)) - (y(18)+y(16)+y(19));
    residual(8) = (y(19)) - (y(20)*y(18));
    residual(9) = (y(11)) - (y(1)^params(5));
    residual(10) = (y(20)) - (x(1));
end
