function real_stage_calc(z_F, x_D, x_B, R_by_Rm, eta)
%Calculates the number of stages required with a given R/Rm ratio for
%distillation of methanol-acetone feed using McCabe Thiele Method

%   z_F is the feed mole fraction
%   x_D is the distillate mole fraction
%   x_B is the bottoms mole fraction
%   R_by_Rm is the value of R/Rm
%   eta is the murphree efficiency in fraction

%To plot equilbrium curve
x_A = 0:0.05:0.8;
y_A = yeq(x_A);
hold on
grid on
axis equal
axis([0 0.8 0 0.8])
plot(x_A, y_A);
plot([0 0.8],[0 0.8]);
xlabel('x_A');
ylabel('y_A');
title('Actual number of stages (\eta = 0.7)');

%Taking specific cases, z_F = 0.35, x_D = 0.75, x_B = 0.001
%D-line
plot([x_D, x_D], [0, x_D], '--');
%B-line
plot([x_B, x_B], [0, x_B], '--');

%Calculations for the recitfying line
l = 0.3;
v = 1 - l;

m_q = -(1 - v)/v;
x_f = z_F;
y_f = x_f;
c_q = y_f - m_q*x_f;
x_0 = 0.25;%initial guess to solve q-line and eq curve
y_fun = @(x_var)(m_q*x_var + c_q - yeq(x_var));
options = optimoptions('fsolve', 'Display', 'off');
x_ep = fsolve(y_fun, x_0, options);%end point of the q-line
y_ep = m_q*x_ep + c_q;
%Feed line
plot([z_F, z_F, x_ep], [0, z_F, y_ep], '--');
Rm = (x_D - y_ep)/(y_ep - x_ep);
R = R_by_Rm*Rm;

%Rectifying line:
m_r = R/(R + 1);
c_r = x_D - m_r*x_D;%intercept of the R line
x_t = 0.25;%temporary end point at R line
y_t = m_r*x_t + c_r;
plot([x_t x_D], [y_t x_D], 'r');

%intersection points for stripping section
x_str = (c_q - c_r)/(m_r - m_q);
y_str = m_r*x_str + c_r;
plot([x_B x_str], [x_B y_str], 'r');
m_str = (y_str - x_B)/(x_str - x_B);
c_str = x_B - m_str*x_B;

%minimum stage

stages = 0;
x_e = x_D;
y_e = x_e;
stage_plot = [];
while (x_e > x_B)
    stage_plot = [stage_plot; [x_e, y_e]];
    y_temp = stage_plot(end, 2);
    x_temp = xeq(y_temp);
    x_e = x_e - eta*(x_e - x_temp);%THIS LINE IS FOR MURPHREE EFF
    stages = stages + 1;
    stage_plot = [stage_plot; [x_e, y_e]];
    y_e = min(m_r*x_e + c_r, m_str*x_e + c_str);
end
stages = stages - (x_e - x_B)/(x_e - stage_plot(end-1, 1));
plot(stage_plot(:, 1), stage_plot(:, 2));
fprintf('Actual number of stages is %.2f \n', stages);
end

function y_A = yeq(x_A)
alpha = -1.5497*x_A + 2.2126;
y_A = (alpha.*x_A)./(1 - x_A + alpha.*x_A);
end

function x_A = xeq(y_A)
alpha = -0.95*(y_A.^2) - 0.8454*y_A + 2.2573;
x_A = y_A./(alpha - alpha.*y_A + y_A);
end