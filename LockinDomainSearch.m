clear all;
%warning off all;
clc;
tau_1=0.5;
tau_2=0.25;
omega_e_free=0;
k_vco=1;

PLLSys = @(t, z, tau_1, tau_2, omega_e_free)([omega_e_free - k_vco * (...
    1 / tau_1 * z(2) + (tau_1 * tau_2 - 1) / tau_1 * z(3) + sin(z(1))...
    ),sin(z(1)),- tau_1 * z(3) + sin(z(1))]');
sys = @(t,z) PLLSys(t, z, tau_1, tau_2, omega_e_free);
Event = @(t, z, omega_e_free,tau_1, k_vco, initial_theta_e)deal([1*((abs(z(1) - round(z(1)/2/pi)*2*pi) ...
    + abs(z(2)-tau_1*omega_e_free/k_vco) + abs(z(3)) >= 1.e-4) && (abs(initial_theta_e - z(1)) >= 2*pi)),1,0]);
len = 200;
view(-2,21);


hold on;
step=0.5;
for x_1_initial=-3.1:step:3
    for x_2_initial=-6:step:6.2
        for theta_e_initial=-pi:step:pi
            min=0;max=0;
            
            xoverFcn = @(t, z) Event(t, z, omega_e_free, tau_1, k_vco, theta_e_initial);
            options = odeset('RelTol', 1.e-3, 'AbsTol', 1.e-3, 'events', xoverFcn);
            [T,Y] = ode15s(sys, [0 len], [theta_e_initial x_1_initial x_2_initial], options);         
            min=theta_e_initial;
            max=theta_e_initial;
            
            for i=1:length(T)
                if(Y(i,1)>max) max=Y(i,1);
                end
                if(Y(i,1)<min) min=Y(i,1);
                end
            end
            
            if(max-min<2*pi)
                scatter3(theta_e_initial, x_1_initial, x_2_initial, 'green', 'LineWidth', 1);
                plot3(Y(:, 1), Y(:, 2), Y(:, 3), 'blue', 'LineWidth', 1);
            end
        end
    end
end

xlabel('\textbf{$\theta_e$}','Interpreter','latex', 'fontsize',20);
ylabel('\textbf{$x_1$}','Interpreter','latex', 'fontsize',20,'rotation',0);
zlabel('\textbf{$x_2$}','Interpreter','latex', 'fontsize',20,'rotation',0);
grid on;
clear T Y;
