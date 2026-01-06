function [dtheta, ds, dx, dy] = f(theta, s, D)
    dtheta = (s^2 - cos(theta)) / s;
    ds = -sin(theta) - D*s^2;
    
    dx = s*cos(theta);
    dy = s*sin(theta);
end

function [theta_new, s_new, x_new, y_new] = rk4_step(theta, s, x, y, D, h)

    [k1_theta, k1_s, k1_x, k1_y] = f(theta, s, D);

    th2 = theta + 0.5*h*k1_theta;
    s2  = s + 0.5*h*k1_s;
    [k2_theta, k2_s, k2_x, k2_y] = f(th2, s2, D);

    th3 = theta + 0.5*h*k2_theta;
    s3  = s + 0.5*h*k2_s;
    [k3_theta, k3_s, k3_x, k3_y] = f(th3, s3, D);

    th4 = theta + h*k3_theta;
    s4  = s + h*k3_s;
    [k4_theta, k4_s, k4_x, k4_y] = f(th4, s4, D);

    theta_new = theta + (h/6)*(k1_theta + 2*k2_theta + 2*k3_theta + k4_theta);
    s_new = s + (h/6)*(k1_s     + 2*k2_s     + 2*k3_s     + k4_s);
    x_new = x + (h/6)*(k1_x     + 2*k2_x     + 2*k3_x     + k4_x);
    y_new = y + (h/6)*(k1_y     + 2*k2_y     + 2*k3_y     + k4_y);
end


function [s_imp, theta_imp, x_imp] = solve_until_impact(theta0, s0, D, h, maxSteps)

    theta = theta0;
    s = s0;
    x = 0;
    y = 1;

    for n = 1:maxSteps
        [theta, s, x, y] = rk4_step(theta, s, x, y, D, h);

        if y <= 0
            s_imp = s;
            theta_imp = theta;
            x_imp = x;
            return
        end
    end

    % Safety fallback (should not happen)
    s_imp = NaN;
    theta_imp = NaN;
    x_imp = NaN;
end

%Accumulated variables re-used for multiples figures
Nth = 80;
theta0_vals = linspace(-pi/3, pi/3, Nth);
D = 0.5;
h = 0.01;
maxSteps = 5000;
vcrit = 0.69;
s0 = 4.5;
x_impact_vals = zeros(size(theta0_vals));
v_y_impact_vals = zeros(size(theta0_vals));

for k = 1:length(theta0_vals)

    theta0 = theta0_vals(k);
    [s_imp, theta_imp, x_imp] = solve_until_impact(theta0, s0, D, h, maxSteps);
    x_impact_vals(k) = x_imp;
    v_y_impact_vals(k) = abs(s_imp * sin(theta_imp));
end

safe = v_y_impact_vals < vcrit;
unsafe = v_y_impact_vals >= vcrit;


figure; hold on;
%xticks([-pi/3 -pi/6 0 pi/6 pi/3])
%xticklabels({'-\pi/3','-\pi/6','0','\pi/6','\pi/3'})

plot(theta0_vals, x_impact_vals, 'Color', [0.7 0.7 0.7], 'Linewidth', 2.5);
x_unsafe = x_impact_vals;
x_unsafe(safe) = NaN;
plot(theta0_vals, x_unsafe, 'r', 'LineWidth', 2);

x_safe = x_impact_vals;
x_safe(unsafe) = NaN;
plot(theta0_vals, x_safe, 'g', 'LineWidth', 2);

xlabel('Launch angle');
ylabel('Horizontal distance and impact');
title('Distance vs Launch angle');
grid on;

legend('Trajectory', 'Egg breaks', 'Egg survives', 'Location', 'best');


figure; hold on;
plot(theta0_vals, v_y_impact_vals, 'LineWidth', 2);
yline(vcrit, '--k', 'LineWidth', 1.5);
%xticks([-pi/3 -pi/6 0 pi/6 pi/3])
%xticklabels({'-\pi/3','-\pi/6','0','\pi/6','\pi/3'})

xlabel('Launch angle \theta_0');
ylabel('Vertical impact speed |v_y|');
title('Impact speed vs launch angle');
grid on;

legend('|v_y^{impact}|', 'Critical speed', 'Location', 'best');
[max(x_impact_vals), max(x_safe)]

fprintf('Unsafe trajectories: %d out of %d\n', sum(unsafe), length(unsafe));
