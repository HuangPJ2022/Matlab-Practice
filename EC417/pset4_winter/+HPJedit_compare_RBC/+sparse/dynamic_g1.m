function [g1, T_order, T] = dynamic_g1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 9
    T_order = -1;
    T = NaN(11, 1);
end
[T_order, T] = HPJedit_compare_RBC.sparse.dynamic_g1_tt(y, x, params, steady_state, T_order, T);
g1_v = NaN(38, 1);
g1_v(1)=(-(exp(x(1))*getPowerDeriv(y(1),params(5),1)));
g1_v(2)=(-(T(5)*y(11)*y(19)*T(10)));
g1_v(3)=(-(1-y(20)));
g1_v(4)=(-(T(5)*y(11)*params(4)*y(19)*T(11)));
g1_v(5)=(-(T(8)*y(11)*(1-params(4))*y(19)*T(10)));
g1_v(6)=(-(T(4)*T(5)));
g1_v(7)=1;
g1_v(8)=(-(T(5)*params(4)*T(7)));
g1_v(9)=(-(T(8)*T(4)*(1-params(4))));
g1_v(10)=1;
g1_v(11)=1;
g1_v(12)=(-(getPowerDeriv(y(13),params(7),1)));
g1_v(13)=1;
g1_v(14)=getPowerDeriv(y(14),T(1),1)/T(1)/(T(2)*T(2));
g1_v(15)=1;
g1_v(16)=(-(y(11)*T(4)*T(9)));
g1_v(17)=(-(y(11)*params(4)*T(7)*T(9)));
g1_v(18)=(-(T(4)*y(11)*(1-params(4))*getPowerDeriv(y(14),(-params(4)),1)));
g1_v(19)=1;
g1_v(20)=(-1);
g1_v(21)=(-1);
g1_v(22)=1;
g1_v(23)=1;
g1_v(24)=(-1)/(T(2)*T(2));
g1_v(25)=(-1);
g1_v(26)=(-(T(5)*y(11)*y(5)*T(10)));
g1_v(27)=(-(T(5)*y(11)*params(4)*y(5)*T(11)));
g1_v(28)=(-(params(2)*params(6)*getPowerDeriv(y(19),params(6)-1,1)));
g1_v(29)=(-(T(8)*y(11)*(1-params(4))*y(5)*T(10)));
g1_v(30)=(-(params(2)*getPowerDeriv(y(19),params(6),1)));
g1_v(31)=y(5);
g1_v(32)=1;
g1_v(33)=(-(params(1)*y(29)/T(3)));
g1_v(34)=(-((-(params(1)*(1+y(22)*y(29)-y(30))*(-(getPowerDeriv(y(24),T(1),1)/T(1)))))/(T(3)*T(3))));
g1_v(35)=(-((-(params(1)*(1+y(22)*y(29)-y(30))))/(T(3)*T(3))));
g1_v(36)=(-(params(1)*y(22)/T(3)));
g1_v(37)=(-((-params(1))/T(3)));
g1_v(38)=(-T(6));
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 10, 31);
end
