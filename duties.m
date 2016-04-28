function duties(Ref)
%Plots the duties as a function of R

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

F = (x_D - x_B)/(z_F - x_B);

Rm = (x_D - y_ep)/(y_ep - x_ep);
reb_duty = [];
cond_duty = [];


for R_by_Rm = Ref    
    
    R = R_by_Rm*Rm;
    V = R + 1;
    V_b = V - F*v;
    lambda = 13860;
    reb_duty = [reb_duty, V_b*lambda];
    R_a = R/((lambda/(28.49*1.8*(55.2 - 30))) + 1);
    ret_ref = R - R_a;
    vap_cont = ret_ref + 1;
    cond_dut = vap_cont*(lambda + 28.49*1.8*(55.2 - 30));
    cond_duty = [cond_duty, cond_dut];
end
hold
plot (Ref, reb_duty, Ref, -cond_duty);
csvwrite('duties.csv',[Ref, reb_duty, -cond_duty]);
grid
xlabel('R/R_m');
ylabel('Duties (Btu/lb mol)');
legend('Reboiler duty', 'Condenser Duty', 'Location', 'East');
%csvwrite('duties.csv', ((cond_duty + reb_duty) - min(cond_duty + reb_duty))/(max(cond_duty + reb_duty) - min(cond_duty + reb_duty)));
end


function y_A = yeq(x_A)
alpha = -1.5497*x_A + 2.2126;
y_A = (alpha.*x_A)./(1 - x_A + alpha.*x_A);
end

function x_A = xeq(y_A)
alpha = -0.95*(y_A.^2) - 0.8454*y_A + 2.2573;
x_A = y_A./(alpha - alpha.*y_A + y_A);
end
