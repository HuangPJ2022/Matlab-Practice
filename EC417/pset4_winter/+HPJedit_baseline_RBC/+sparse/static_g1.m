function [g1, T_order, T] = static_g1(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 8
    T_order = -1;
    T = NaN(9, 1);
end
[T_order, T] = HPJedit_baseline_RBC.sparse.static_g1_tt(y, x, params, T_order, T);
g1_v = NaN(33, 1);
g1_v(1)=(-(T(2)*T(3)));
g1_v(2)=1-exp(x(1))*getPowerDeriv(y(1),params(5),1);
g1_v(3)=(-(T(3)*params(4)*T(4)));
g1_v(4)=(-(T(5)*T(2)*(1-params(4))));
g1_v(5)=(-(T(1)*params(1)*y(9)));
g1_v(6)=1;
g1_v(7)=1;
g1_v(8)=(-(T(1)*T(6)));
g1_v(9)=1;
g1_v(10)=1;
g1_v(11)=(-(y(1)*T(2)*T(7)));
g1_v(12)=(-(y(1)*params(4)*T(4)*T(7)));
g1_v(13)=(-(T(2)*y(1)*(1-params(4))*getPowerDeriv(y(4),(-params(4)),1)));
g1_v(14)=(-(T(3)*y(1)*y(9)*T(8)));
g1_v(15)=1-(1-y(10));
g1_v(16)=(-(T(3)*y(1)*params(4)*y(9)*T(9)));
g1_v(17)=(-(T(5)*y(1)*(1-params(4))*y(9)*T(8)));
g1_v(18)=(-1);
g1_v(19)=(-1);
g1_v(20)=1;
g1_v(21)=1;
g1_v(22)=(-1)/(y(8)*y(8))-params(1)*(1+y(2)*y(9)-y(10))*(-1)/(y(8)*y(8));
g1_v(23)=(-(T(6)*(-y(3))/(y(8)*y(8))));
g1_v(24)=(-1);
g1_v(25)=(-(T(1)*params(1)*y(2)));
g1_v(26)=(-(T(3)*y(1)*y(5)*T(8)));
g1_v(27)=(-(T(3)*y(1)*params(4)*y(5)*T(9)));
g1_v(28)=(-(params(2)*params(6)*getPowerDeriv(y(9),params(6)-1,1)));
g1_v(29)=(-(T(5)*y(1)*(1-params(4))*y(5)*T(8)));
g1_v(30)=(-(params(2)*getPowerDeriv(y(9),params(6),1)));
g1_v(31)=(-(T(1)*(-params(1))));
g1_v(32)=y(5);
g1_v(33)=1;
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 10, 10);
end
