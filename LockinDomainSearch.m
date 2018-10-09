clear all;
%warning off all;
clc;
tau_1=0.5;
tau_2=0.25;
omega_e_free=0;

sys = @(t,z) PLLSys(t, z, tau_1, tau_2, omega_e_free);
len = 20;
%xoverFcn = @(t, z) Event(t, z);
%options = odeset('RelTol', 1.e-3, 'AbsTol', 1.e-3, 'events', xoverFcn);
view(-2,21);


hold on;
step=1;
for x_1_initial=-3.1:step:3
    for x_2_initial=-6:step:6.2
        for theta_e_initial=-pi:step:pi
            min=0;max=0;
            xoverFcn = @(t, z) Event(t, z, theta_e_initial);
            options = odeset('RelTol', 1.e-3, 'AbsTol', 1.e-3, 'events', xoverFcn);
            [T,Y] = ode15s(sys, [0 len], [theta_e_initial x_1_initial x_2_initial], options);
            %min_y = min(Y(:,1,1));
            %max_y = max(Y(:,1,1));
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

function [value,isterminal,direction] = Event( t, z, initial_theta_e)

% 1 - theta_e
% 2 - x_1
% 3 - x_2
theta_e = z(1);
x1 = z(2);
x2 = z(3);
is_in_equilibria = abs(theta_e - round(theta_e/2/pi)*2*pi) ...
    + abs(x1) + abs(x2) < 1.e-4;
is_cycle_slipping = abs(initial_theta_e - theta_e) > 2*pi;
if (is_in_equilibria || is_cycle_slipping)
    value = 0;
else 
    value = 1;
end
isterminal=1;
direction=0;

end

function dz = PLLSys(t, z,tau_1,tau_2,omega_e_free)

dz = zeros(3,1);

% 1 - theta_e
% 2 - x_1
% 3 - x_2
theta_e = z(1);
x1 = z(2);
x2 = z(3);


%TODO pass this from arguments
k_vco = 1;

dx1 = sin(theta_e);
dx2 = - tau_1 * x2 + sin(theta_e);
dtheta_e = omega_e_free - k_vco * (...
    1 / tau_1 * x1 + (tau_1 * tau_2 - 1) / tau_1 * x2 + sin(theta_e)...
    );

dz(1) = dtheta_e;
dz(2) = dx1;
dz(3) = dx2;

end

