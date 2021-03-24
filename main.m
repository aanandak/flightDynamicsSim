% clc
% clear all
% close all

% Find trim conditions with given initial conditions
V = 250;
h = 10000;

[theta0, delta_e, T] = findTrim(V, h);


%% Evaluate stability of trim conditions
stability = stabAnalysis(V, h, theta0, delta_e, T);


%% Simulate flight dynamics

U0 = V*cos(theta0);
W0 = V*sin(theta0);


%s = [x  z  u   w   theta   q  eta1 eta2]  format
s0 = [0; h; U0; W0; theta0; 0; 0; 0]; 
t = 0:0.1:5000;

% TECS PID
error_T = 0;
error_L = 0;
count = 0;
h_dot_old = 0;
v_dot_old = 0;

%; error_T; error_L; h_dot_old; v_dot_old
f = [delta_e; T];
[tt, xx1] = ode45(@(t,x)FW_longitudinal_dynamics(t, x, f), t, s0);


%% Plotting
subplot(3,2,1)
plot(tt, xx1(:, 1));
title('x-dist');

subplot(3,2,2)
plot(tt, xx1(:, 2));
title('Altitude');

subplot(3,2,3)
plot(tt, xx1(:, 3));
title('Horiz Speed');

subplot(3,2,4)
plot(tt, xx1(:, 4));
title('Vert Speed');

subplot(3,2,5)
plot(tt, xx1(:, 5)*180/pi);
title('Theta');

subplot(3,2,6)
plot(tt, xx1(:, 6));
title('q');

for pp = 1:6
    subplot(3,2,pp)
    hold on
    grid on
end
