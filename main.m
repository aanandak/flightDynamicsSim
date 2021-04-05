% clc
% clear all
% close all

% Find trim conditions with given initial conditions
V = 250;
h = 10000;

[theta0, delta_e, T] = findTrim(V, h);


% Trim at desired conditions

V2 = 260;
h2 = 12000;

[theta02, delta_e2, T2] = findTrim(V2, h2);

%% Evaluate stability of trim conditions
% stability = stabAnalysis(V, h, theta0, delta_e, T);


%% Simulate flight dynamics

U0 = V*cos(theta0);
W0 = V*sin(theta0);

h_sp = 15000;
v_sp = 280;


% TECS PID
error_T = 0;
error_L = 0;


% format: s =[x  z  u  w  theta  q ] 
s0 = [0; h; U0; W0; theta0; 0; error_T; error_L]; 
t = 0:2500;


f = [delta_e; T; h_sp; v_sp];
[tt, xx1] = ode45(@(t,x)FW_longitudinal_dynamics(t, x, f), t, s0);


%% Plotting
subplot(4,2,1)
plot(tt, xx1(:, 1));
title('x-dist');

subplot(4,2,2)
plot(tt, xx1(:, 2));
title('Altitude');

subplot(4,2,3)
plot(tt, xx1(:, 3));
title('Horiz Speed');

subplot(4,2,4)
plot(tt, xx1(:, 4));
title('Vert Speed');

subplot(4,2,5)
plot(tt, xx1(:, 5)*180/pi);
title('Theta');

subplot(4,2,6)
plot(tt, xx1(:, 6));
title('q');

subplot(4,2,7)
plot(tt, xx1(:, 7));
title('Thrust Integral');

subplot(4,2,8)
plot(tt, xx1(:, 8));
title('Elevator Integral');

for pp = 1:8
    subplot(4,2,pp)
    hold on
    grid on
end


final_v = sqrt(xx1(end, 3)^2 + xx1(end, 4)^2)
final_h = xx1(end, 2)
