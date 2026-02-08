function [g1, T_order, T] = static_g1(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T_order, T)
if nargin < 8
    T_order = -1;
    T = NaN(3, 1);
end
[T_order, T] = HPJedit_baseline_RBC.sparse.static_g1_tt(y, x, params, T_order, T);
g1_v = NaN(28, 1);
g1_v(1)=(-(T(2)*T(3)));
g1_v(2)=1-getPowerDeriv(y(1),params(5),1);
g1_v(3)=(-T(1));
g1_v(4)=1;
g1_v(5)=(-(1/((1+y(10))*y(8))));
g1_v(6)=1;
g1_v(7)=params(6)/((1-y(4))*(1-y(4)));
g1_v(8)=(-(y(1)*T(2)*getPowerDeriv(y(4),1-params(4),1)));
g1_v(9)=(-((-(y(7)*(1-params(4))))/(y(4)*y(4))));
g1_v(10)=(-(T(3)*y(1)*getPowerDeriv(y(5),params(4),1)));
g1_v(11)=(-((-(y(7)*params(4)))/(y(5)*y(5))));
g1_v(12)=1-(1-params(2));
g1_v(13)=(-1);
g1_v(14)=(-1);
g1_v(15)=1;
g1_v(16)=(-((1-params(4))/y(4)));
g1_v(17)=(-(params(4)/y(5)));
g1_v(18)=1;
g1_v(19)=(-(1+y(10)))/((1+y(10))*y(8)*(1+y(10))*y(8))-(1+y(2)-params(2))*1/(1+params(1))*(-(1+y(10)))/((1+y(10))*y(8)*(1+y(10))*y(8));
g1_v(20)=(-((-((1+y(10))*y(3)))/((1+y(10))*y(8)*(1+y(10))*y(8))));
g1_v(21)=(-1);
g1_v(22)=(-y(10));
g1_v(23)=(-1);
g1_v(24)=1;
g1_v(25)=(-y(8))/((1+y(10))*y(8)*(1+y(10))*y(8))-(1+y(2)-params(2))*1/(1+params(1))*(-y(8))/((1+y(10))*y(8)*(1+y(10))*y(8));
g1_v(26)=(-((-(y(8)*y(3)))/((1+y(10))*y(8)*(1+y(10))*y(8))));
g1_v(27)=(-y(8));
g1_v(28)=1;
if ~isoctave && matlab_ver_less_than('9.8')
    sparse_rowval = double(sparse_rowval);
    sparse_colval = double(sparse_colval);
end
g1 = sparse(sparse_rowval, sparse_colval, g1_v, 10, 10);
end
