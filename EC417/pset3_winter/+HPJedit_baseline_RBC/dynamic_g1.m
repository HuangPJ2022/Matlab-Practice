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
g1 = zeros(10, 16);
g1(1,13)=(-T(1));
g1(1,10)=(-(1+y(12)))/((1+y(12))*y(10)*(1+y(12))*y(10));
g1(1,14)=(-((1+y(13)-params(2))*1/(1+params(1))*(-(1+y(15)))/((1+y(15))*y(14)*(1+y(15))*y(14))));
g1(1,12)=(-y(10))/((1+y(12))*y(10)*(1+y(12))*y(10));
g1(1,15)=(-((1+y(13)-params(2))*1/(1+params(1))*(-y(14))/((1+y(15))*y(14)*(1+y(15))*y(14))));
g1(2,5)=(-(1/((1+y(12))*y(10))));
g1(2,6)=params(6)/((1-y(6))*(1-y(6)));
g1(2,10)=(-((-((1+y(12))*y(5)))/((1+y(12))*y(10)*(1+y(12))*y(10))));
g1(2,12)=(-((-(y(10)*y(5)))/((1+y(12))*y(10)*(1+y(12))*y(10))));
g1(3,3)=(-(T(2)*T(3)));
g1(3,6)=(-(y(3)*T(2)*getPowerDeriv(y(6),1-params(4),1)));
g1(3,2)=(-(T(3)*y(3)*getPowerDeriv(y(2),params(4),1)));
g1(3,9)=1;
g1(4,5)=1;
g1(4,6)=(-((-(y(9)*(1-params(4))))/(y(6)*y(6))));
g1(4,9)=(-((1-params(4))/y(6)));
g1(5,4)=1;
g1(5,2)=(-((-(y(9)*params(4)))/(y(2)*y(2))));
g1(5,9)=(-(params(4)/y(2)));
g1(6,2)=(-(1-params(2)));
g1(6,7)=1;
g1(6,8)=(-1);
g1(7,8)=(-1);
g1(7,9)=1;
g1(7,10)=(-1);
g1(7,11)=(-1);
g1(8,10)=(-y(12));
g1(8,11)=1;
g1(8,12)=(-y(10));
g1(9,1)=(-(getPowerDeriv(y(1),params(5),1)));
g1(9,3)=1;
g1(10,12)=1;
g1(10,16)=(-1);

end
