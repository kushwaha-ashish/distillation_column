function min_stages(z_F, x_D, x_B)
%This function calculates the minimum number of stages required for
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
plot([0 0.8],[0 0.8]);
xlabel('x_A');
ylabel('y_A');
title('number of stages for acetone-methanol distillation');

%Taking specific cases, z_F = 0.35, x_D = 0.75, x_B = 0.001

%D-line
plot([x_D, x_D], [0, x_D], '--');
%B-line
plot([x_B, x_B], [0, x_B], '--');
%F-line
plot([z_F, z_F], [0, z_F], '--');

%minimum stage
stages = 0;
x_e = x_D;
stage_plot = [];
while (x_e > x_B)
    stage_plot = [stage_plot; [x_e, x_e]];
    y_temp = stage_plot(end, 2);
    x_temp = xeq(y_temp);
    stages = stages + 1;
    stage_plot = [stage_plot; [x_temp, y_temp]];
    x_e = x_temp;
end
stages = stages - (x_e - x_B)/(x_e - stage_plot(end-1, 1));
plot(stage_plot(:, 1), stage_plot(:, 2));
fprintf('Minimum number of stages is %.2f \n', stages);
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