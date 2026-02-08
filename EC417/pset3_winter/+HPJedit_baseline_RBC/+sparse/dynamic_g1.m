function [g1, T_order, T] = dynamic_g1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 9
    T_order = -1;
    T = NaN(3, 1);
end
[T_order, T] = HPJedit_baseline_RBC.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
g1_v = NaN(33, 1);
g1_v(1)=(-(getPowerDeriv(y(1),params(5),1)));
g1_v(2)=(-(T(3)*y(11)*getPowerDeriv(y(5),params(4),1)));
g1_v(3)=(-((-(y(17)*params(4)))/(y(5)*y(5))));
g1_v(4)=(-(1-params(2)));
g1_v(5)=(-(T(2)*T(3)));
g1_v(6)=1;
g1_v(7)=1;
g1_v(8)=(-(1/((1+y(20))*y(18))));
g1_v(9)=1;
g1_v(10)=params(6)/((1-y(14))*(1-y(14)));
g1_v(11)=(-(y(11)*T(2)*getPowerDeriv(y(14),1-params(4),1)));
g1_v(12)=(-((-(y(17)*(1-params(4))))/(y(14)*y(14))));
g1_v(13)=1;
g1_v(14)=(-1);
g1_v(15)=(-1);
g1_v(16)=1;
g1_v(17)=(-((1-params(4))/y(14)));
g1_v(18)=(-(params(4)/y(5)));
g1_v(19)=1;
g1_v(20)=(-(1+y(20)))/((1+y(20))*y(18)*(1+y(20))*y(18));
g1_v(21)=(-((-((1+y(20))*y(13)))/((1+y(20))*y(18)*(1+y(20))*y(18))));
g1_v(22)=(-1);
g1_v(23)=(-y(20));
g1_v(24)=(-1);
g1_v(25)=1;
g1_v(26)=(-y(18))/((1+y(20))*y(18)*(1+y(20))*y(18));
g1_v(27)=(-((-(y(18)*y(13)))/((1+y(20))*y(18)*(1+y(20))*y(18))));
g1_v(28)=(-y(18));
g1_v(29)=1;
g1_v(30)=(-T(1));
g1_v(31)=(-((1+y(22)-params(2))*1/(1+params(1))*(-(1+y(30)))/((1+y(30))*y(28)*(1+y(30))*y(28))));
g1_v(32)=(-((1+y(22)-params(2))*1/(1+params(1))*(-y(28))/((1+y(30))*y(28)*(1+y(30))*y(28))));
g1_v(33)=(-1);
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 10, 31);
end
