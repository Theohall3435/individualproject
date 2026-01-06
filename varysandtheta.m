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
    y = 1.0;

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


Ns = 40;
Nth = 40;
s0_vals = linspace(0.15, 2.2, Ns);
theta0_vals = linspace(-pi/3, pi/3, Nth);
D = 0.5; %will vary later
h = 0.01;
maxSteps = 5000;
vcrit = 0.69;
survival = zeros(Nth, Ns);

for i = 1:Nth
    for j = 1:Ns

        theta0 = theta0_vals(i);
        s0 = s0_vals(j);

        [s_imp, theta_imp, ~] = solve_until_impact(theta0, s0, D, h, maxSteps);

        v_y_imp = abs(s_imp * sin(theta_imp));

        survival(i,j) = (v_y_imp < vcrit);
    end
end

figure;
imagesc(s0_vals, theta0_vals, survival);
set(gca, 'YDir', 'normal');
colormap([1 0 0; 0 1 0]);
caxis([0 1]);
xlabel('Initial speed s_0');
ylabel('Initial angle \theta_0');
title('Egg survival region for fixed D');
