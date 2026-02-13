function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = HPJedit_baseline_RBC.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 9
    T = [T; NaN(9 - size(T, 1), 1)];
end
T(6) = getPowerDeriv(y(3)/y(8),params(7),1);
T(7) = getPowerDeriv(y(4),1-params(4),1);
T(8) = getPowerDeriv(y(9)*y(5),params(4),1);
T(9) = getPowerDeriv(y(9)*y(5),params(4)-1,1);
end
