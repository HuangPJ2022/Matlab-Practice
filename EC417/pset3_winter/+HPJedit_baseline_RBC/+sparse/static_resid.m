function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(3, 1);
end
[T_order, T] = HPJedit_baseline_RBC.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(10, 1);
    residual(1) = (1/((1+y(10))*y(8))) - (T(1)*(1+y(2)-params(2)));
    residual(2) = (params(6)/(1-y(4))) - (y(3)/((1+y(10))*y(8)));
    residual(3) = (y(7)) - (y(1)*T(2)*T(3));
    residual(4) = (y(3)) - (y(7)*(1-params(4))/y(4));
    residual(5) = (y(2)) - (y(7)*params(4)/y(5));
    residual(6) = (y(5)) - (y(5)*(1-params(2))+y(6));
    residual(7) = (y(7)) - (y(8)+y(6)+y(9));
    residual(8) = (y(9)) - (y(10)*y(8));
    residual(9) = (y(1)) - (y(1)^params(5));
    residual(10) = (y(10)) - (x(1));
end
