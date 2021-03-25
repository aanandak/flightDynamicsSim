function sdot = FW_longitudinal_dynamics(~, s, f) 
    % s = [x z u w theta q]
    x = s(1);
    z = s(2);
    u = s(3);
    w = s(4);
    theta = s(5);
    q = s(6);

    Data;
    qp = 0.5*rho*(u^2+w^2);

    alpha = atan2(w, u); 
    ih = alpha; 
    delta_ec = f(1);
    gamma = theta - alpha;
    g = 9.81;
    
   
    % TECS params
    v_sp = 250;
    h_sp = 10000;
    gamma_sp = 0;
    gamma_tilda = gamma - gamma_sp;

    v = sqrt(u^2+w^2); 
    V_tilda = (v - v_sp);
    h = z;
    
    eta1 = s(7);
    eta2 = s(8);
    
    % TECS PID Constants
    KV = 0.3;
    KP = -0.005;
    KD = 1;
    
    KTP = 0.3;
    KTI = 0.05;
  
    KEP = 0.025;
    KEI = 0.005;


    % Eq. 111
    d_thr_bar = (KTP/(KTP+1))*(KV*V_tilda/g + sin(gamma_tilda) + sin(gamma))...
        + (KTI/(KTP+1))*eta1;
    
    % Eq. 112
    d_elv_bar = KP*theta + KD*q +...
        KEP*(KV*V_tilda/g - sin(gamma_tilda) + sin(gamma) - d_thr_bar)...
        + KEI*eta2;

    
    delta_ec = delta_ec + d_elv_bar;
    delta_thr = f(2) + d_thr_bar;

    
    % Dynamics Equations
    cL = Cl0 + Cla*alpha + Clih*ih + Clde*delta_ec;
    cDa = 2*cL*Cla/(pi *ARw*e);
    cD = Cd_0 + cDa*alpha;
    cM = Cm0 + Cma*alpha + Cmih*ih + Cmde*delta_ec;
    
    L = qp * cL * Sw;
    D = qp * cD * Sw;
    Moment = cM * qp * Sw * Cw;
    X = -D*cos(alpha) + L*sin(alpha) + delta_thr;
    Z = -D*sin(alpha) - L*cos(alpha);

    xdot = u*cos(theta) + w*sin(theta);
    zdot = -u*sin(theta) + w*cos(theta);
    udot = -q*w - g*sin(theta) + X/m;
    wdot = q*u + g*cos(theta) + Z/m; 
    thetadot = q;
    qdot = Moment/ Iyy;
    
    
    % Eq 113
    eta1_dot = KV*V_tilda/g + sin(gamma_tilda) + sin(gamma) - d_thr_bar;

    % Eq 114
    eta2_dot = KV*V_tilda/g - sin(gamma_tilda) + sin(gamma) - d_thr_bar;
    
    sdot = [xdot; zdot; udot; wdot; thetadot; qdot; eta1_dot; eta2_dot];
end
