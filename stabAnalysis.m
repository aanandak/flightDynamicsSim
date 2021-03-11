function stability = stabAnalysis(V, h, theta0, delta_ec, delta_thr)
    syms z u w t q;

    
    % Lift, drag, moment coefs
    Data;
    % Conditions
    sigma   = lambda^4.2561;   

    
   

    Cw = 7.49;                 
    Thr_maxsl = 4*138800;      
    g = 9.81;               
    Iyy = 30513547.0;       
    m = GW/g;               


    alpha = t; 
    ih = t;

    qp      = 0.5*rho*(u^2+w^2);   
    

    % Dynamics Equations
    cL = Cl0 + Cla*alpha + Clih*ih + Clde*delta_ec;  % Lift Coef
    cDa = 2*cL*Cla/(pi *ARw*e);                 
    cD = Cd_0 + cDa*alpha;                      % Drag Coef
    cM = Cm0 + Cma*alpha + Cmih*ih + Cmde*delta_ec;  % Moment Coef


    L = qp * cL * Sw;      % Lift 
    D = qp * cD * Sw;      % Drag
%     Thr = Thr_maxsl * sigma * delta_thr;      % Thrust
    Moment = cM * qp * Sw * Cw;         % Moment
%     X = -D*cos(alpha) - L*sin(alpha) + Thr*cos(t);
    X = -D*cos(alpha) + L*sin(alpha) + delta_thr;
    Z = -D*sin(alpha) - L*cos(alpha);% - Thr*sin(t); % This extra term made no difference


    zdot = -u*sin(t) + w*cos(t);
    udot = -q*w - g*sin(t) + X/m;
    wdot = q*u + g*cos(t) + Z/m; 
    tdot = q;
    qdot = Moment/ Iyy;


    A = [ diff(zdot, z) diff(zdot, u) diff(zdot, w) diff(zdot, t) diff(zdot, q);
          diff(udot, z) diff(udot, u) diff(udot, w) diff(udot, t) diff(udot, q);
          diff(wdot, z) diff(wdot, u) diff(wdot, w) diff(wdot, t) diff(wdot, q);
          diff(tdot, z) diff(tdot, u) diff(tdot, w) diff(tdot, t) diff(tdot, q);
          diff(qdot, z) diff(qdot, u) diff(qdot, w) diff(qdot, t) diff(qdot, q)];


    A = subs(A, q, 0);
    A = subs(A, z, h);
    A = subs(A, t, theta0);
    A = subs(A, u, V*cos(theta0));
    A = subs(A, w, V*sin(theta0));


    stability = vpa(eig(A), 3);

end

