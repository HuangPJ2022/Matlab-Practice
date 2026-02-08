function [y, T, residual, g1] = static_3(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(8, 1);
  T(1)=y(1)*y(5)^params(4);
  T(2)=y(4)^(1-params(4));
  residual(1)=(y(7))-(T(1)*T(2));
  residual(2)=(y(3))-(y(7)*(1-params(4))/y(4));
  residual(3)=(y(2))-(y(7)*params(4)/y(5));
  residual(4)=(y(5))-(y(5)*(1-params(2))+y(6));
  residual(5)=(y(7))-(y(8)+y(6)+y(9));
  residual(6)=(y(9))-(y(10)*y(8));
  residual(7)=(1/((1+y(10))*y(8)))-(1/((1+y(10))*y(8))*1/(1+params(1))*(1+y(2)-params(2)));
  residual(8)=(params(6)/(1-y(4)))-(y(3)/((1+y(10))*y(8)));
if nargout > 3
    g1_v = NaN(22, 1);
g1_v(1)=(-(T(1)*getPowerDeriv(y(4),1-params(4),1)));
g1_v(2)=(-((-(y(7)*(1-params(4))))/(y(4)*y(4))));
g1_v(3)=params(6)/((1-y(4))*(1-y(4)));
g1_v(4)=1;
g1_v(5)=(-(1/((1+y(10))*y(8))));
g1_v(6)=(-(T(2)*y(1)*getPowerDeriv(y(5),params(4),1)));
g1_v(7)=(-((-(y(7)*params(4)))/(y(5)*y(5))));
g1_v(8)=1-(1-params(2));
g1_v(9)=(-1);
g1_v(10)=(-1);
g1_v(11)=1;
g1_v(12)=(-((1-params(4))/y(4)));
g1_v(13)=(-(params(4)/y(5)));
g1_v(14)=1;
g1_v(15)=(-1);
g1_v(16)=1;
g1_v(17)=1;
g1_v(18)=(-(1/((1+y(10))*y(8))*1/(1+params(1))));
g1_v(19)=(-1);
g1_v(20)=(-y(10));
g1_v(21)=(-(1+y(10)))/((1+y(10))*y(8)*(1+y(10))*y(8))-(1+y(2)-params(2))*1/(1+params(1))*(-(1+y(10)))/((1+y(10))*y(8)*(1+y(10))*y(8));
g1_v(22)=(-((-((1+y(10))*y(3)))/((1+y(10))*y(8)*(1+y(10))*y(8))));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 8, 8);
end
end
