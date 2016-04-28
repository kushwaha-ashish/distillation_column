function min_reflux(z_F, x_D, x_B)
%This function calculates the minimum reflux required for
%distillation of methanol-acetone feed using McCabe Thiele Method

%   z_F is the feed mole fraction
%   x_D is the distillate mole fraction
%   x_B is the bottoms mole fraction

%To plot equilbrium curve
x_A = 0:0.05:0.8;
y_A = yeq(x_A);
hold on
grid on
axis equal
axis([0 0.8 0 0.8])
plot(x_A, y_A);
plot([0 0.8],[0 0.8], '--');%x=y line
xlabel('x_A');
ylabel('y_A');
title('Pinch point and rectifying line for acetone-methanol distillation');

%Taking specific cases, z_F = 0.35, x_D = 0.75, x_B = 0.001

%D-line
plot([x_D, x_D], [0, x_D], '--');
%B-line
plot([x_B, x_B], [0, x_B], '--');

%Feed line
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
y_ep = m_q*x_ep + c_q;%end point of the q-line
plot([z_F, z_F, x_ep], [0, z_F, y_ep], '--');
plot([x_D x_ep], [x_D y_ep], 'r');
% minimum reflux from the pinch point
r = (x_D - y_ep)/(y_ep - x_ep);
fprintf('Pinch point is x = %.4f, y = %.4f \n', x_ep, y_ep);
fprintf('The minimum reflux required is %.2f \n', r);
end

%Function to calculate y_eq for a given x
function y_A = yeq(x_A)
alpha = -1.5497*x_A + 2.2126;
y_A = (alpha.*x_A)./(1 - x_A + alpha.*x_A);
end

%Function to calculate x_eq for a given y
function x_A = xeq(y_A)
alpha = -0.95*(y_A.^2) - 0.8454*y_A + 2.2573;
x_A = y_A./(alpha - alpha.*y_A + y_A);
end