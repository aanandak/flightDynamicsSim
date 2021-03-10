%% Find trim conditions with given initial conditions
V = 250;
h = 10000;

[theta0, delta_e, delta_thr] = findTrim(V, h);


%% Evaluate stability of trim conditions
stability = stabAnalysis(V, h, theta0, delta_e, delta_thr);


%% Simulate flight dynamics

% theta = 6.177426627134937e-02;
U0 = V*cos(theta0);
W0 = V*sin(theta0);


% format: s =[x  z  u  w  theta  q ] 
s0 = [0; h; U0; W0; theta0; 0;]; 
t = 0:0.1:10000;
u_in  = 0.7;
[tt, xx1] = ode45(@(t,x)FW_longitudinal_dynamics(t, x), t, s0);


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