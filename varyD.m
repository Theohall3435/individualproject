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
    y = 20.0;

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
s0_vals = linspace(0.2, 3.0, Ns);
theta0_vals = linspace(-pi/3, pi/3, Nth);
D = 0.5;
h = 0.01;
maxSteps = 5000;
vcrit = 0.69;
survival = zeros(Nth, Ns);
theta = -pi/2 + 0.05;
s0 = 4.5;
D_values = linspace(0, 10, 50);
v_y_impact_values = zeros(size(D_values));

for k = 1:length(D_values)

    D = D_values(k);

    [s_imp, theta_imp, ~] = solve_until_impact(theta, s0, D, h, maxSteps);

    v_y_impact_values(k) = abs(s_imp * sin(theta_imp));
end

for i = 1:Nth
    for j = 1:Ns

        theta0 = theta0_vals(i);
        s0 = s0_vals(j);
        D_map = 0.5;
        [s_imp, theta_imp, ~] = solve_until_impact(theta0, s0, D, h, maxSteps);

        v_y_imp = abs(s_imp * sin(theta_imp));

        survival(i,j) = (v_y_imp < vcrit);
    end
end

safe = v_y_impact_values < vcrit;
unsafe = v_y_impact_values >= vcrit;


idx = find(v_y_impact_values < vcrit, 1, 'first');
if ~isempty(idx)
    D_star = D_values(idx);
    fprintf('Critical drag D* ≈ %.2f\n', D_star);
end

figure; hold on;
plot(D_values, v_y_impact_values, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
%unsafe segment
v_unsafe = v_y_impact_values;
v_unsafe(safe) = NaN;
plot(D_values, v_unsafe, 'r', 'LineWidth', 2);

%safe segment
v_safe = v_y_impact_values;
v_safe(unsafe) = NaN;
plot(D_values, v_safe, 'g', 'LineWidth', 2);

% Critical threshold
yline(vcrit, '--k', 'LineWidth', 1.5);

xlabel('Drag parameter D');
ylabel('Vertical impact speed |v_y|');
title('Egg survival vs drag parameter');
grid on;

fprintf('max |v_y| = %.3f\n', max(v_y_impact_values));
