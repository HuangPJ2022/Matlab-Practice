function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = HPJedit_baseline_RBC.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 10
    T = [T; NaN(10 - size(T, 1), 1)];
end
T(7) = getPowerDeriv(y(13)/y(18),params(7),1);
T(8) = getPowerDeriv(y(14),1-params(4),1);
T(9) = getPowerDeriv(y(19)*y(5),params(4),1);
T(10) = getPowerDeriv(y(19)*y(5),params(4)-1,1);
end
