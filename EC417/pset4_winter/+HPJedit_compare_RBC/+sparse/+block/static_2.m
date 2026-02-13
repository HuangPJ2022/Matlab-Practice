function [y, T, residual, g1] = static_2(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(9, 1);
  residual(1)=(y(4))-(y(3)^params(7));
  T(2)=(y(9)*y(5))^params(4);
  T(3)=y(4)^(1-params(4));
  residual(2)=(y(7))-(y(1)*T(2)*T(3));
  residual(3)=(y(7))-(y(8)+y(6));
  residual(4)=(y(5))-(y(6)+y(5)*(1-y(10)));
  T(4)=y(8)-y(4)^(1+1/params(7))/(1+1/params(7));
  residual(5)=(1/T(4))-(params(1)*(1+y(2)*y(9)-y(10))/T(4));
  T(5)=y(1)*params(4)*(y(9)*y(5))^(params(4)-1);
  residual(6)=(y(2))-(T(3)*T(5));
  residual(7)=(y(2))-(params(2)*params(6)*y(9)^(params(6)-1));
  T(6)=y(4)^(-params(4));
  residual(8)=(y(3))-(T(2)*y(1)*(1-params(4))*T(6));
  residual(9)=(y(10))-(params(2)*y(9)^params(6));
  T(7)=getPowerDeriv(y(4),1+1/params(7),1)/(1+1/params(7));
  T(8)=getPowerDeriv(y(4),1-params(4),1);
  T(9)=getPowerDeriv(y(9)*y(5),params(4),1);
  T(10)=getPowerDeriv(y(9)*y(5),params(4)-1,1);
if nargout > 3
    g1_v = NaN(29, 1);
g1_v(1)=(-(getPowerDeriv(y(3),params(7),1)));
g1_v(2)=1;
g1_v(3)=1;
g1_v(4)=1;
g1_v(5)=(-1);
g1_v(6)=(-1)/(T(4)*T(4))-(-(params(1)*(1+y(2)*y(9)-y(10))))/(T(4)*T(4));
g1_v(7)=(-1);
g1_v(8)=(-1);
g1_v(9)=1;
g1_v(10)=(-(y(1)*T(2)*T(8)));
g1_v(11)=T(7)/(T(4)*T(4))-(-(params(1)*(1+y(2)*y(9)-y(10))*(-T(7))))/(T(4)*T(4));
g1_v(12)=(-(T(5)*T(8)));
g1_v(13)=(-(T(2)*y(1)*(1-params(4))*getPowerDeriv(y(4),(-params(4)),1)));
g1_v(14)=(-(T(3)*y(1)*y(5)*T(9)));
g1_v(15)=(-(params(1)*y(2)/T(4)));
g1_v(16)=(-(T(3)*y(1)*params(4)*y(5)*T(10)));
g1_v(17)=(-(params(2)*params(6)*getPowerDeriv(y(9),params(6)-1,1)));
g1_v(18)=(-(T(6)*y(1)*(1-params(4))*y(5)*T(9)));
g1_v(19)=(-(params(2)*getPowerDeriv(y(9),params(6),1)));
g1_v(20)=(-(params(1)*y(9)/T(4)));
g1_v(21)=1;
g1_v(22)=1;
g1_v(23)=(-(T(3)*y(1)*y(9)*T(9)));
g1_v(24)=1-(1-y(10));
g1_v(25)=(-(T(3)*y(1)*params(4)*y(9)*T(10)));
g1_v(26)=(-(T(6)*y(1)*(1-params(4))*y(9)*T(9)));
g1_v(27)=y(5);
g1_v(28)=(-((-params(1))/T(4)));
g1_v(29)=1;
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 9, 9);
end
end
