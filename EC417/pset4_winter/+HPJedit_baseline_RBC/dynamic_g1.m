function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = HPJedit_baseline_RBC.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(10, 17);
g1(1,13)=(-(T(1)*params(1)*y(15)));
g1(1,10)=(-1)/(y(10)*y(10));
g1(1,14)=(-(params(1)*(1+y(13)*y(15)-y(16))*(-1)/(y(14)*y(14))));
g1(1,15)=(-(T(1)*params(1)*y(13)));
g1(1,16)=(-(T(1)*(-params(1))));
g1(2,5)=(-(1/y(10)*T(7)));
g1(2,6)=1;
g1(2,10)=(-(T(7)*(-y(5))/(y(10)*y(10))));
g1(3,3)=(-(T(2)*T(3)));
g1(3,6)=(-(y(3)*T(2)*T(8)));
g1(3,2)=(-(T(3)*y(3)*y(11)*T(9)));
g1(3,9)=1;
g1(3,11)=(-(T(3)*y(3)*y(2)*T(9)));
g1(4,8)=(-1);
g1(4,9)=1;
g1(4,10)=(-1);
g1(5,2)=(-(1-y(12)));
g1(5,7)=1;
g1(5,8)=(-1);
g1(5,12)=y(2);
g1(6,1)=(-(exp(x(it_, 1))*getPowerDeriv(y(1),params(5),1)));
g1(6,3)=1;
g1(6,17)=(-T(4));
g1(7,3)=(-(T(3)*params(4)*T(5)));
g1(7,4)=1;
g1(7,6)=(-(y(3)*params(4)*T(5)*T(8)));
g1(7,2)=(-(T(3)*y(3)*params(4)*y(11)*T(10)));
g1(7,11)=(-(T(3)*y(3)*params(4)*y(2)*T(10)));
g1(8,4)=1;
g1(8,11)=(-(params(2)*params(6)*getPowerDeriv(y(11),params(6)-1,1)));
g1(9,3)=(-(T(6)*T(2)*(1-params(4))));
g1(9,5)=1;
g1(9,6)=(-(T(2)*y(3)*(1-params(4))*getPowerDeriv(y(6),(-params(4)),1)));
g1(9,2)=(-(T(6)*y(3)*(1-params(4))*y(11)*T(9)));
g1(9,11)=(-(T(6)*y(3)*(1-params(4))*y(2)*T(9)));
g1(10,11)=(-(params(2)*getPowerDeriv(y(11),params(6),1)));
g1(10,12)=1;

end
