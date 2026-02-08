function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(11)=y(1)^params(5);
  y(20)=x(1);
end
