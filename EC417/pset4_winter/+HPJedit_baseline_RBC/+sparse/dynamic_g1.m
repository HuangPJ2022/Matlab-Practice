function [g1, T_order, T] = dynamic_g1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 9
    T_order = -1;
    T = NaN(10, 1);
end
[T_order, T] = HPJedit_baseline_RBC.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
g1_v = NaN(37, 1);
g1_v(1)=(-(exp(x(1))*getPowerDeriv(y(1),params(5),1)));
g1_v(2)=(-(T(3)*y(11)*y(19)*T(9)));
g1_v(3)=(-(1-y(20)));
g1_v(4)=(-(T(3)*y(11)*params(4)*y(19)*T(10)));
g1_v(5)=(-(T(6)*y(11)*(1-params(4))*y(19)*T(9)));
g1_v(6)=(-(T(2)*T(3)));
g1_v(7)=1;
g1_v(8)=(-(T(3)*params(4)*T(5)));
g1_v(9)=(-(T(6)*T(2)*(1-params(4))));
g1_v(10)=1;
g1_v(11)=1;
g1_v(12)=(-(1/y(18)*T(7)));
g1_v(13)=1;
g1_v(14)=1;
g1_v(15)=(-(y(11)*T(2)*T(8)));
g1_v(16)=(-(y(11)*params(4)*T(5)*T(8)));
g1_v(17)=(-(T(2)*y(11)*(1-params(4))*getPowerDeriv(y(14),(-params(4)),1)));
g1_v(18)=1;
g1_v(19)=(-1);
g1_v(20)=(-1);
g1_v(21)=1;
g1_v(22)=1;
g1_v(23)=(-1)/(y(18)*y(18));
g1_v(24)=(-(T(7)*(-y(13))/(y(18)*y(18))));
g1_v(25)=(-1);
g1_v(26)=(-(T(3)*y(11)*y(5)*T(9)));
g1_v(27)=(-(T(3)*y(11)*params(4)*y(5)*T(10)));
g1_v(28)=(-(params(2)*params(6)*getPowerDeriv(y(19),params(6)-1,1)));
g1_v(29)=(-(T(6)*y(11)*(1-params(4))*y(5)*T(9)));
g1_v(30)=(-(params(2)*getPowerDeriv(y(19),params(6),1)));
g1_v(31)=y(5);
g1_v(32)=1;
g1_v(33)=(-(T(1)*params(1)*y(29)));
g1_v(34)=(-(params(1)*(1+y(22)*y(29)-y(30))*(-1)/(y(28)*y(28))));
g1_v(35)=(-(T(1)*params(1)*y(22)));
g1_v(36)=(-(T(1)*(-params(1))));
g1_v(37)=(-T(4));
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 10, 31);
end
