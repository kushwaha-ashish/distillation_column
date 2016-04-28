function reflux_effect(Ref)
%Plotting the effect of reflux on the number of stages
eta = 0.7;
z_F = 0.35;
x_D = 0.75;
x_B = 0.001;

l = 0.3;
v = 1 - l;

m_q = -(1 - v)/v;
x_f = z_F;
y_f = x_f;
c_q = y_f - m_q*x_f;
x_0 = 0.25;%initial guess to solve q-line and eq curve
y_fun = @(x_var)(m_q*x_var + c_q - yeq(x_var));
x_ep = fsolve(y_fun, x_0);
y_ep = m_q*x_ep + c_q;

Rm = (x_D - y_ep)/(y_ep - x_ep);
stage_arr = [];

for R_by_Rm = Ref
    R = R_by_Rm*Rm;
    
    %Rectifying line:
    m_r = R/(R + 1);
    c_r = x_D - m_r*x_D;%intercept of the R line
    x_t = 0.25;%temporary end point at R line
    y_t = m_r*x_t + c_r;

    %intersection points for stripping section
    x_str = (c_q - c_r)/(m_r - m_q);
    y_str = m_r*x_str + c_r;
    m_str = (y_str - x_B)/(x_str - x_B);
    c_str = x_B - m_str*x_B;

    stages = 0;
    x_e = x_D;
    y_e = x_e;
    stage_plot = [];
    
    while (x_e > x_B)
        stage_plot = [stage_plot; [x_e, y_e]];
        y_temp = stage_plot(end, 2);
        x_temp = xeq(y_temp);
        x_temp = x_e - eta*(x_e - x_temp);%THIS LINE IS FOR MURPHREE EFF
        stage_plot = [stage_plot; [x_temp, y_temp]];
        stages = stages + 1;
        x_e = x_temp;
        y_e = min(m_r*x_e + c_r, m_str*x_e + c_str);
    end
    
    stages = stages - (x_e - x_B)/(x_e - stage_plot(end-1, 1));

    plot(stage_plot(:, 1), stage_plot(:, 2));


    stage_arr = [stage_arr stages];
end
plot(Ref, stage_arr);
xlabel('R/R_m');
ylabel('Number of stages');
grid on
end

function y_A = yeq(x_A)
alpha = -1.5497*x_A + 2.2126;
y_A = (alpha.*x_A)./(1 - x_A + alpha.*x_A);
end

function x_A = xeq(y_A)
alpha = -0.95*(y_A.^2) - 0.8454*y_A + 2.2573;
x_A = y_A./(alpha - alpha.*y_A + y_A);
end

