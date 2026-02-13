function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = HPJedit_baseline_RBC.static_resid_tt(T, y, x, params);
end
residual = zeros(10, 1);
    residual(1) = (T(1)) - (T(1)*params(1)*(1+y(2)*y(9)-y(10)));
    residual(2) = (y(4)) - ((y(3)/y(8))^params(7));
    residual(3) = (y(7)) - (y(1)*T(2)*T(3));
    residual(4) = (y(7)) - (y(8)+y(6));
    residual(5) = (y(5)) - (y(6)+y(5)*(1-y(10)));
    residual(6) = (y(1)) - (y(1)^params(5)*exp(x(1)));
    residual(7) = (y(2)) - (T(3)*y(1)*params(4)*T(4));
    residual(8) = (y(2)) - (params(2)*params(6)*y(9)^(params(6)-1));
    residual(9) = (y(3)) - (T(2)*y(1)*(1-params(4))*T(5));
    residual(10) = (y(10)) - (params(2)*y(9)^params(6));

end
