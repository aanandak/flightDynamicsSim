function stability = stabAnalysis(V, h, theta0, delta_ec, delta_thr)
    syms z u w t q;

    % Conditions
    rhosl   = 1.225;              
    lambda  = (1-22.556e-6*z); % Using this changes makes first eig non-zero
    % lambda  = 0.7744;   
    rho     = rhosl*lambda^4.2561;
    qp      = 0.5*rho*(u^2+w^2);   
    sigma   = lambda^4.2561;   

    Sw = 363.12;               
    ARw = 10.03;               
    Cw = 7.49;                 
    Thr_maxsl = 4*138800;      
    g = 9.81;               
    GW = 2500000;           
    Iyy = 30513547.0;       
    m = GW/g;               

    % Lift, drag, moment coefs
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


    alpha = t; 
    ih = t;


    % Dynamics Equations
    cL = Cl0 + Cla*alpha + Clih*ih + Clde*delta_ec;  % Lift Coef
    cDa = 2*cL*Cla/(pi *ARw*e);                 
    cD = Cd_0 + cDa*alpha;                      % Drag Coef
    cM = Cm0 + Cma*alpha + Cmih*ih + Cmde*delta_ec;  % Moment Coef


    L = qp * cL * Sw;      % Lift 
    D = qp * cD * Sw;      % Drag
    Thr = Thr_maxsl * sigma * delta_thr;      % Thrust
    Moment = cM * qp * Sw * Cw;         % Moment
    X = -D*cos(alpha) - L*sin(alpha) + Thr*cos(t);
    Z = D*sin(alpha) - L*cos(alpha);% - Thr*sin(t); % This extra term made no difference


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

