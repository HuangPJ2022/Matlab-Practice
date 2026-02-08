function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(6, 1);
  y(13)=y(17)*(1-params(4))/y(14);
  y(19)=y(20)*y(18);
  T(1)=y(11)*y(5)^params(4);
  T(2)=y(14)^(1-params(4));
  residual(1)=(y(17))-(T(1)*T(2));
  residual(2)=(y(17))-(y(18)+y(16)+y(19));
  residual(3)=(params(6)/(1-y(14)))-(y(13)/((1+y(20))*y(18)));
  residual(4)=(y(15))-(y(5)*(1-params(2))+y(16));
  residual(5)=(y(12))-(y(17)*params(4)/y(5));
  residual(6)=(1/((1+y(20))*y(18)))-(1/(1+params(1))*1/((1+y(30))*y(28))*(1+y(22)-params(2)));
if nargout > 3
    g1_v = NaN(18, 1);
g1_v(1)=(-(T(2)*y(11)*getPowerDeriv(y(5),params(4),1)));
g1_v(2)=(-(1-params(2)));
g1_v(3)=(-((-(y(17)*params(4)))/(y(5)*y(5))));
g1_v(4)=1;
g1_v(5)=1;
g1_v(6)=(-((1-params(4))/y(14)/((1+y(20))*y(18))));
g1_v(7)=(-(params(4)/y(5)));
g1_v(8)=(-1);
g1_v(9)=(-1);
g1_v(10)=(-(T(1)*getPowerDeriv(y(14),1-params(4),1)));
g1_v(11)=params(6)/((1-y(14))*(1-y(14)))-(-(y(17)*(1-params(4))))/(y(14)*y(14))/((1+y(20))*y(18));
g1_v(12)=1;
g1_v(13)=1;
g1_v(14)=(-(1+y(20)));
g1_v(15)=(-((-((1+y(20))*y(13)))/((1+y(20))*y(18)*(1+y(20))*y(18))));
g1_v(16)=(-(1+y(20)))/((1+y(20))*y(18)*(1+y(20))*y(18));
g1_v(17)=(-(1/(1+params(1))*1/((1+y(30))*y(28))));
g1_v(18)=(-((1+y(22)-params(2))*1/(1+params(1))*(-(1+y(30)))/((1+y(30))*y(28)*(1+y(30))*y(28))));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 6, 18);
end
end
