function sdot = FW_longitudinal_dynamics(~, s, f) 
    % s = [x z u w theta q]
    x = s(1);
    z = s(2);
    u = s(3);
    w = s(4);
    theta = s(5);
    q = s(6);

    Data;
    qp = 0.5*rho*(u^2+w^2);           % dynamic pressure

    alpha = theta; 
    ih = theta; 
    delta_ec = f(1);
    

    % Dynamics Equations
    cL = Cl0 + Cla*alpha + Clih*ih + Clde*delta_ec;
    cDa = 2*cL*Cla/(pi *ARw*e);
    cD = Cd_0 + cDa*alpha;
    cM = Cm0 + Cma*alpha + Cmih*ih + Cmde*delta_ec;
    
    L = qp * cL * Sw;
    D = qp * cD * Sw;
    Moment = cM * qp * Sw * Cw;
    X = -D*cos(alpha) + L*sin(alpha) + f(2);
    Z = -D*sin(alpha) - L*cos(alpha);

    xdot = u*cos(theta) + w*sin(theta);
    zdot = -u*sin(theta) + w*cos(theta);
    udot = -q*w - g*sin(theta) + X/m;
    wdot = q*u + g*cos(theta) + Z/m; 
    thetadot = q;
    qdot = Moment/ Iyy;
    
    
    sdot = [xdot; zdot; udot; wdot; thetadot; qdot];
end

