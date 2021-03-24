function [theta, delta_e, T] = findTrim(V, h)
Data;
%% Eq 68

t0 = [0, pi/3];
% myfun = @(t) (cos(t)*(Cl0+(Cla+Clih)*t+Clde*(-(Cm0+(Cma+Cmih)*t)/Cmde)) + sin(t)*(Cd_0+2*(Cl0+Cla*t)*Cla/(pi *ARw*e)*t) - GW*cos(t)/(C1*V^2));
myfun = @(t) (cos(t)*(Cl0+(Cla+Clih)*t+Clde*(-(Cm0+(Cma+Cmih)*t)/Cmde)) + sin(t)*(Cd_0 + ((2*(Cl0+(Cla+Clih)*t+Clde*(-(Cm0+(Cma+Cmih)*t)/Cmde))*Cla)/(pi*ARw*e))*t) - GW*cos(t)/(C1*V^2));

x = fzero(myfun, t0);

theta = x;
delta_e = -(Cm0 + (Cma+Cmih)*x)/Cmde;
alpha = theta;
ih = alpha;

cL = Cl0 + Cla*alpha + Clih*ih + Clde*delta_e;
cDa = 2*cL*Cla/(pi *ARw*e);
cD = Cd_0 + cDa*alpha;
T = (C1 * V^2 * cD * sec(theta));


end
