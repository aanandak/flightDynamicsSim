function sdot = FW_longitudinal_dynamics(~, s) 
    % s = [x z u w theta q]
    x = s(1);
    z = s(2);
    u = s(3);
    w = s(4);
    theta = s(5);
    q = s(6);

    % Conditions
    rhosl   = 1.225;               % Air density sea level
%     lambda  = (1-22.556e-6*z);
    lambda = 0.7744;
    rho     = rhosl*lambda^4.2561;  % actual density
    qp      = 0.5*rho*(u^2+w^2);           % dynamic pressure
    sigma   = lambda^4.2561;       % air density coef
    Sw = 363.12;                % wing area
    ARw = 10.03;                % wing aspect ratio
    Cw = 7.49;                  % wing chord length
    Thr_maxsl = 4*138800;       % max Thrust at sea level 
    g = 9.81;               % Grav accel
    GW = 2500000;           % Gross Weight, N
    Iyy = 30513547.0;       % kg-m2
    m = GW/g;               % Mass, kg
    
    % Lift, drag, moment coef components
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
    
%     theta = 0.06177;
%     alpha0 = 3.7355 * pi/180;
    alpha = theta; %alpha0*0 + atan2(w,u); 
    ih = theta; %-4.642563651*pi/180*0 + theta*0;
    delta_ec = -4.927003209569403e-01; 
    delta_thr = 2689;
    
    % Dynamics Equations
    cL = Cl0 + Cla*alpha + Clih*ih + Clde*delta_ec;  % Lift Coef
    cDa = 2*cL*Cla/(pi *ARw*e);                 
    cD = Cd_0 + cDa*alpha;                      % Drag Coef
    cM = Cm0 + Cma*alpha + Cmih*ih + Cmde*delta_ec;  % Moment Coef
    
    L = qp * cL * Sw;      % Lift 
    D = qp * cD * Sw;      % Drag
    Thr = Thr_maxsl * sigma * delta_thr;      % Thrust
    Moment = cM * qp * Sw * Cw;         % Moment
    X = -D*cos(alpha) - L*sin(alpha) + Thr*cos(alpha);
    Z = D*sin(alpha) - L*cos(alpha) - Thr*sin(theta); % positive downwards

    
    xdot = u*cos(theta) + w*sin(theta);
    zdot = -u*sin(theta) + w*cos(theta);
    udot = -q*w - g*sin(theta) + X/m;
    wdot = q*u + g*cos(theta) + Z/m; 
    thetadot = q;
    qdot = Moment/ Iyy;
    
   
    
    sdot = [xdot; zdot; udot; wdot; thetadot; qdot];
end