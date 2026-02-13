function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  T(1)=exp(x(1));
  y(11)=y(1)^params(5)*T(1);
end
