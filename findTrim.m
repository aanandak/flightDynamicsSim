function [theta, delta_e, delta_thr] = findTrim(V, h)

% System Details
Cla = 5.9598;
Cl0 = 0.2301;
Clih = 0.8299;
Clde = 0.2391;
e = 0.85;
Cm0 = -0.0812;
Cma = -3.1069;
Cmih = -3.4077;
Cmde = -0.9816;
Cd_0 = 0.0172;

Sw = 363.12;
ARw = 10.03; 
GW = 2500000; 
rhosl = 1.225;              
lambda = (1-22.556e-6*h);
rho = rhosl*lambda^4.2561;
C1 = 0.5*rho*Sw;




%% Eq 68

t0 = [0, pi/3];
myfun = @(t) (cos(t)*(Cl0+(Cla+Clih)*t+Clde*(-(Cm0+(Cma+Cmih)*t)/Cmde)) + sin(t)*(Cd_0+2*(Cl0+Cla*t)*Cla/(pi *ARw*e)*t) - GW*cos(t)/(C1*V^2));
x = fzero(myfun, t0);


theta = x;
delta_e = -(Cm0 + (Cma+Cmih)*x)/Cmde;
alpha = theta;
ih = theta;

cL = Cl0 + Cla*alpha + Clih*ih + Clde*delta_e;  % Lift Coef
cDa = 2*cL*Cla/(pi *ARw*e);                 
cD = Cd_0 + cDa*alpha;                      % Drag Coef
T = (C1 * V^2 * cD * sec(theta)); 
Thr_maxsl = 4*138800;  

delta_thr = T/Thr_maxsl;

end
