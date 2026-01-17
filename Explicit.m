clear
clc
close all

U = 1;                
nu = 1.4e-5; 
den = 1.225;
L = 8;                 
H = 0.06;             
dx = 0.001;            
dy = 0.001;          
Nx = int64(L/dx + 1);         
Ny = int64(H/dy + 1);         

u = zeros(Ny, Nx);     
v = zeros(Ny, Nx);    

u(:, 1) = U;           
u(1, :) = 0;           
v(1, :) = 0;           
v(end, :) = 0;        
u(end, :) = U;         

x = linspace(0, L, Nx);
y = linspace(0, H, Ny);


for i = 1:Nx-1
    fprintf('%d %d\n', i, Nx);

    for j = 2:Ny-1
        u(j, i+1) = u(j, i) + nu * (u(j+1, i) - 2 * u(j, i) + u(j-1, i)) / u(j, i) ...
                    * dx / (dy^2) - ((u(j+1, i) - u(j-1, i)) / 2 )* (v(j, i) / u(j, i)) * (dx / dy);
        
        v(j, i+1) = v(j-1, i+1)  - 0.5 * ((u(j,i+1) - u(j,i)) / dx + (u(j-1,i+1) - u(j-1,i)) / dx) * dy;
    end
end

% Blasius analytical solutions
Re = U * x / nu;
delta_blausius = 5 * x ./ sqrt(Re);
delta_blausius(1)=0;
delta_star_blausius = 1.72 * x ./ sqrt(Re);
delta_star_blausius(1)=0;
delta_2_blausius = 0.664 * x ./ sqrt(Re);
delta_2_blausius(1)=0;
cf_blausius = 0.664 ./ sqrt(Re);
%cf_blausius(1)=0;

% Numerical solutions
[~, idx] = max(u >= 0.99 * U, [], 1); 
delta_numerical = y(idx);

f = 1 - u / U;
delta_star_numerical = trapz(y, f);

f2 = (u / U) .* (1 - u / U);
delta_2_numerical = trapz(y, f2);

t_wall = nu * (u(2, :) / dy);
cf_numerical = t_wall / (0.5 * den * U^2);
%%
%PLOTTING STUFF
figure
plot(x, delta_blausius, 'r-', 'LineWidth', 1.5)
hold on
plot(x, delta_numerical, 'b--', 'LineWidth', 1.5)
xlabel('x (m)')
ylabel('Boundary Layer Thickness (m)')
title('Boundary Layer Thickness Comparison')
legend('Blasius Solution', 'Numerical Solution')
grid on
hold off

figure
plot(x, delta_star_blausius, 'r-', 'LineWidth', 1.5)
hold on
plot(x, delta_star_numerical, 'b--', 'LineWidth', 1.5);
xlabel('x (m)')
ylabel('Boundary Layer Displacement Thickness (m)')
title('Boundary Layer Displacement Thickness Comparison')
legend('Blasius Solution', 'Numerical Solution')
grid on
hold off

figure
plot(x, delta_2_blausius, 'r-', 'LineWidth', 1.5); hold on
plot(x, delta_2_numerical, 'b--', 'LineWidth', 1.5)
xlabel('x (m)')
ylabel('Boundary Layer Momentum Loss Thickness (m)')
title('Boundary Layer Momentum Loss Comparison')
legend('Blasius Solution', 'Numerical Solution')
grid on
hold off

figure
plot(x, cf_blausius, 'r-', 'LineWidth', 1.5); hold on
plot(x, cf_numerical, 'b--', 'LineWidth', 1.5)
xlabel('x (m)')
ylabel('Friction Coefficient')
title('Friction Coefficient Comparison')
legend('Blasius Solution', 'Numerical Solution')
grid on
hold off


% Error calculations
% Boundary Layer Thickness
error_delta_mae = mean(abs(delta_blausius - delta_numerical));
error_delta_rmse = sqrt(mean((delta_blausius - delta_numerical).^2));
error_delta_max = max(abs(delta_blausius - delta_numerical));

% Displacement Thickness
error_delta_star_mae = mean(abs(delta_star_blausius - delta_star_numerical));
error_delta_star_rmse = sqrt(mean((delta_star_blausius - delta_star_numerical).^2));
error_delta_star_max = max(abs(delta_star_blausius - delta_star_numerical));

% Momentum Thickness
error_delta_2_mae = mean(abs(delta_2_blausius - delta_2_numerical));
error_delta_2_rmse = sqrt(mean((delta_2_blausius - delta_2_numerical).^2));
error_delta_2_max = max(abs(delta_2_blausius - delta_2_numerical));

% Friction Coefficient
error_cf_mae = mean(abs(cf_blausius(2:end) - cf_numerical(2:end)));
error_cf_rmse = sqrt(mean((cf_blausius(2:end) - cf_numerical(2:end)).^2));
error_cf_max = max(abs(cf_blausius(2:end) - cf_numerical(2:end)));

disp('Errors for Boundary Layer Thickness (delta):')
disp(['MAE: ', num2str(error_delta_mae)])
disp(['RMSE: ', num2str(error_delta_rmse)])
disp(['Max Error: ', num2str(error_delta_max)])

disp('Errors for Displacement Thickness (delta*):')
disp(['MAE: ', num2str(error_delta_star_mae)])
disp(['RMSE: ', num2str(error_delta_star_rmse)])
disp(['Max Error: ', num2str(error_delta_star_max)])

disp('Errors for Momentum Thickness (delta2):')
disp(['MAE: ', num2str(error_delta_2_mae)])
disp(['RMSE: ', num2str(error_delta_2_rmse)])
disp(['Max Error: ', num2str(error_delta_2_max)])

disp('Errors for Friction Coefficient (cf):')
disp(['MAE: ', num2str(error_cf_mae)])
disp(['RMSE: ', num2str(error_cf_rmse)])
disp(['Max Error: ', num2str(error_cf_max)])





plot_contours = input(['Contour και  διανύσματα της ταχύτητας eίναι βαριά. Plot? (yes/no): '], 's');

if strcmpi(plot_contours, 'yes')
    [X, Y] = meshgrid(x, y);
    figure
    contourf(X, Y, sqrt(u.^2 + v.^2), 50, 'LineColor', 'none') 
    hold on
    plot(X(:), Y(:), 'k.', 'MarkerSize', 0.5); 
    plot(x, delta_blausius, 'r-', 'LineWidth', 3)
    plot(x, delta_numerical, 'b--', 'LineWidth', 3)
    colorbar;
    xlabel('x (m)')
    ylabel('y (m)')
    title('Velocity |V| Contour Plot')
    grid on
    hold off

    figure
    quiver(X, Y, u, v, 'k')
    yline(0, 'b', 'LineWidth', 0.5)
    xlabel('x (m)')
    ylabel('y (m)')
    title('Velocity Vector Field')
    axis equal
    grid on
else
    disp('Plots were skipped by user request.');
end





