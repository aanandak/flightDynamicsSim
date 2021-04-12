function sdot = FW_longitudinal_dynamics(~, s, f) 
    % s = [x z u w theta q]
    x = s(1);
    z = s(2);
    u = s(3);
    w = s(4);
    theta = s(5);
    q = s(6);

    Data;
    qp = 0.5*rho*(u^2+w^2);% dynamic pressure

    alpha = atan2(w, u);
    ih = theta; 
    delta_ec = f(1);
    delta_thr = f(2);
    gamma = theta - alpha;
    

    % TECS params
    h_sp = f(3);
    v_sp = f(4);
   
    v = sqrt(u^2+w^2);
    h = z;
    
    gamma_sp = 0.0001*(h_sp - h);
   
    c(1) = m;
    c(2) = v_sp;
    c(3) = h_sp;
    c(4) = gamma_sp;
    c(5) = v;
    c(6) = h;
    c(7) = gamma;
    c(8) = s(7);
    c(9) = s(8);
    
    [d_thr, d_elv, error_T_dot, error_L_dot] = TECS_Controller(c);
    
    delta_ec = delta_ec + d_elv;
    delta_thr = delta_thr + d_thr;
       
    
    % Dynamics Equations
    cL = Cl0 + Cla*alpha + Clih*ih + Clde*delta_ec;
    cDa = 2*cL*Cla/(pi *ARw*e);
    cD = Cd_0 + cDa*alpha;
    cM = Cm0 + Cma*alpha + Cmih*ih + Cmde*delta_ec;
    
    L = qp * cL * Sw;
    D = qp * cD * Sw;
    Moment = cM * qp * Sw * Cw;
    X = -D*cos(alpha) + L*sin(alpha) + delta_thr;%f(2);
    Z = -D*sin(alpha) - L*cos(alpha);

    xdot = u*cos(theta) + w*sin(theta);
    zdot = -u*sin(theta) + w*cos(theta);
    udot = -q*w - g*sin(theta) + X/m;
    wdot = q*u + g*cos(theta) + Z/m; 
    thetadot = q;
    qdot = Moment/ Iyy;
    
    
    sdot = [xdot; zdot; udot; wdot; thetadot; qdot; error_T_dot; error_L_dot];
end




