function phi = incomp_flow(dx, dy)

    % Constants
    w = 1.8;  % relaxation factor
    tol = 1e-5;
    K = 2*(1/dx^2 + 1/dy^2);
    K_airfoil = 2/dx^2 + 1/dy^2;
    max_res = inf;
    max_iter = 20000;
    iter = 0;

    % Domain
    L = 5;
    W = 5;
    x = -L:dx:L;
    y = 0:dy:W;

    % Grid initialization
    phi = zeros(length(x), length(y));
    phi_new = zeros(size(phi));

    while (max_res > tol) && (iter < max_iter)
        iter = iter + 1;

        for i = 2:length(x)-1
            for j = 1:length(y)-1

                % Bottom boundary
                if j == 1
                    if (x(i) > 0) && (x(i) <= 1)
                        f_prime = 0.6 * (0.14845*x(i)^(-0.5) - 0.126 - 0.7032*x(i) + 0.8529*x(i)^2 - 0.406*x(i)^3);
                        phi_new(i, j) = 1/K_airfoil * ((phi(i+1, j) + phi(i-1, j))/dx^2 + ...
                                         (phi(i, j+1)/dy - f_prime)/dy);
                    else
                        phi_new(i, j) = 1/K * ((phi(i+1, j) + phi(i-1, j))/dx^2 + 2*phi(i, j+1)/dy^2);
                    end
                else
                    % SOR update
                    phi_temp = 1/K * ((phi(i+1, j) + phi_new(i-1, j))/dx^2 + ...
                                      (phi(i, j+1) + phi_new(i, j-1))/dy^2);
                    phi_new(i, j) = phi(i, j) + w * (phi_temp - phi(i, j));
                end
            end
        end

        % Compute max residual
        max_res = compute_residual(phi_new, dx, dy);

        % Copy new to old
        phi = phi_new;
    end

    % Final status
    if iter == max_iter
        fprintf('Reached max iterations with max residual = %.3e\n', max_res);
    else
        fprintf('Converged at iteration %d\n', iter);
    end

    % Contour plot
    figure;
    [X, Y] = meshgrid(x, y);
    contourf(X, Y, transpose(phi));
    colorbar;
    xlabel('x');
    ylabel('y');
    title('Potential flow')

end


function max_res = compute_residual(phi, dx, dy)
    res = zeros(size(phi));

    for j = 2:size(phi, 2)-1
        for i = 2:size(phi, 1)-1
            res(i, j) = (phi(i+1, j) - 2*phi(i, j) + phi(i-1, j))/dx^2 + ...
                        (phi(i, j+1) - 2*phi(i, j) + phi(i, j-1))/dy^2;
        end
    end

    max_res = max(abs(res(:)));
end


function surface_vel = compute_surface_vel(phi, dx)

    % Domain
    L = 5;
    x = -L:dx:L;
    x_airfoil = 0:dx:1;
    
    % Initialize solution vector
    surface_vel = zeros(size(x_airfoil));

    % Surface velocity computation
    for i = 2:length(x)-1
        if (x(i) > 0) && (x(i) <= 1)
            surface_vel(round(x(i)/dx)+1) = 1 + (phi(i, 1) - phi(i-1, 1))/dx;
        end
    end
     
    % Experimental data
    x_exp = [0, 0.5, 1.25, 2.5, 5, 7.5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 95, 100] ./ 100;
    surface_vel_exp = [0, 0.64, 1.01, 1.241, 1.378, 1.402, 1.411, 1.411, 1.399, 1.378, 1.350, 1.288, 1.228, 1.166, 1.109, 1.044, 0.956, 0.906, 0];
    
    % Plotting and comparing surface velocity with experimental data
    figure;
    hold on
    plot(x_airfoil, surface_vel.^2);
    plot(x_exp, surface_vel_exp);
    hold off
    legend('SOR', 'Experimental')
    xlabel('x/c')
    ylabel('(v/V)^2')
    title('Surface velocity')
end


%% Incompressible flow and surface velocity computation
clc, clear, close all

dx = 0.02;
dy = 0.02;
phi = incomp_flow(dx, dy);
surface_vel = compute_surface_vel(phi, dx);