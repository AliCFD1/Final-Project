clc
clear
close all;

% User data
n_x = 102;        % number of (physical) cells along x
n_y = n_x;        % number of (physical) cells along y 
vel_in = 1;       % inlet velocity
Re = 500;         % Reynolds number
CFL = 0.5;        % CFL number
L = 1;            % length [m]
L_x = L;          % length [m]
L_y = L;        % length [m]
nu = vel_in*L/Re; % kinematic viscosity [m2/s]
final_time = 2;  % total time of simulation [s]

% Boundary conditions
u_n  = 0; % north wall velocity [m/s]
u_s  = 0; % south wall velocity [m/s]
v_e  = 0; % east wall velocity [m/s]
v_w  = 0; % west wall velocity [m/s]
u_in = 1; % inlet velocity on the west side [m/s]

% Parameters for solver
solver_option  = 'SOR'; % Options: SOR | GMRES | Direct
max_iter       = 50000; % maximum number of iterations
omega          = 1.8;   % SOR coefficient
maximum_error  = 1e-4;  % error for convergence

% Inlet section: west side
n_in_start = 1/2*(n_y+2)  ;   % first cell index 
n_in_end   = 1  *(n_y+2)-1;   % last cell index

% Outlet section: east side
n_out_start = 1           ;  % first cell index 
n_out_end   = 1/2*(n_y+2)+1; % last cell index

% Data processing
if (mod(n_x,2)~=0 || mod(n_y,2)~=0)
    error('Only even number of cells can be accepted (for graphical purposes only)');
end
if (mod(n_y+2,4)~=0)
    error('The total number of cells along y (i.e. ny+2) must be divisible by 4');
end

% Process the grid
h_x = L_x / n_x ; % grid step (uniform grid) [m]
h_y = L_y / n_y ; % grid step (uniform grid) [m]


% Inlet - Outlet section areas
A_in  = h_y * (n_in_end  - n_in_start  + 1); % inlet  section area [m]
A_out = h_y * (n_out_end - n_out_start + 1); % outlet section area [m]

% Estimated max velocity
u_max = max([abs(u_n) , abs(u_s) , abs(v_e) , abs(v_w) , u_in , u_in*A_in/A_out ]); % maximum velocity [m/s]

% Time step
stabilizer = 0.50;                                    % safety factor for time step (stability)
dt_1       = min(h_x,h_y)^2/(4*nu);                   % time step (diffusion stability)  [s]
dt_2       = (4*nu)/u_max^2;                          % time step (convection stability) [s]
dt_CFL     = (CFL*min(h_x,h_y))/vel_in;               % time step (based on CFL number)  [s]
dt         = stabilizer*min(min(dt_1, dt_2),dt_CFL);  % time step (stability)            [s]

n_steps = final_time/dt;   % number of steps

fprintf('Time step size: %.4f\n',  dt);
fprintf('Reynolds number: %.0f\n', Re);

% Grid construction
x     = 0:h_x:L_x;       % grid coordinates (x axis)
y     = 0:h_y:L_y;       % grid coordinates (y axis)
[X,Y] = meshgrid(x,y);   % MATLAB grid

% pre - definition
% Main fields (velocities and pressure)
u_0 = zeros(n_x+1,n_y+2);
v_0 = zeros(n_x+2,n_y+1);
p_0 = zeros(n_x+2,n_y+2);

% Temporary velocity fields
u_t = zeros(n_x+1,n_y+2);
v_t = zeros(n_x+2,n_y+1);

% Fields used only for graphical post-processing purposes
u_u = zeros(n_x+1,n_y+1);
v_v = zeros(n_x+1,n_y+1);
p_p = zeros(n_x+1,n_y+1);

% Coefficient for pressure equation
alpha              = zeros(n_x+2,n_y+2)+h_x*h_y/(2*h_x^2+2*h_y^2); % internal cells
alpha(2,3:n_y)     = h_x*h_y/(2*h_x^2+h_y^2);                      % west cells
alpha(n_x+1,3:n_y) = h_x*h_y/(2*h_x^2+h_y^2);                      % east cells
alpha(3:n_x,2)     = h_x*h_y/(h_x^2+2*h_y^2);                      % south cells
alpha(3:n_x,n_y+1) = h_x*h_y/(h_x^2+2*h_y^2);                      % north cells
alpha(2,2)         = h_x*h_y/(h_x^2+h_y^2);                        % corner cells
alpha(2,n_y+1)     = h_x*h_y/(h_x^2+h_y^2);                        % corner cells     
alpha(n_x+1,2)     = h_x*h_y/(h_x^2+h_y^2);                        % corner cells
alpha(n_x+1,n_y+1) = h_x*h_y/(h_x^2+h_y^2);                        % corner cells

% Adjusting alpha coefficients, no correction in prescribed inlet sections, treated like wall cells.
alpha(n_x+1,n_out_start:n_out_end) = h_x*h_y/(2*h_x^2+h_y^2);

% Adapting gamma coefficients for outlets; correction is possible due to unknown outlet velocity. Outlet cells exhibit internal cell behavior
alpha(n_x+1,n_out_start:n_out_end) = h_x*h_y/(2*h_x^2+2*h_y^2);

% Initial conditions
u_0(1,n_in_start:n_in_end)       = u_in; % inlet section: fixed velocity [m/s]
u_0(n_x+1,n_out_start:n_out_end) = u_in; % outlet section: fixed velocity [m/s]
u_0(2:n_x,2:n_y+1)               = u_in; % internal points: fixed velocity [m/s]
u_t                              = u_0;  % temporary velocity [m/s]

% Solution over time
% Initialize figure for plotting Q_out, Q_in, and the difference
figure(11);
xlabel('Iteration');
ylabel('Flow Rate [m^2/s]');
title('Inlet and Outlet Flow Rates vs Iteration');
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);

% Initialize arrays to store Q_in, Q_out, and their difference values at each iteration
Q_in_values   = zeros(1, round(n_steps));
Q_out_values  = zeros(1, round(n_steps));
Q_diff_values = zeros(1, round(n_steps));

t=0.0;
for is=1:n_steps
    
    % Boundary conditions
    u_0(1:n_x+1,1)     = 2*u_s-u_0(1:n_x+1,2);     % south wall
    u_0(1:n_x+1,n_y+2) = 2*u_n-u_0(1:n_x+1,n_y+1); % north wall
    v_0(1,1:n_y+1)     = 2*v_w-v_0(2,1:n_y+1);     % west wall
    v_0(n_x+2,1:n_y+1) = 2*v_e-v_0(n_x+1,1:n_y+1); % east wall
    
    % Over-writing inlet conditions    
    u_0(1,n_in_start:n_in_end) = u_in;               % fixed velocity [m/s]
    
    % Over-writing outlet conditions
    u_0(n_x+1,n_out_start:n_out_end) = u_0(n_x,n_out_start:n_out_end);    % zero-gradient     
    v_0(n_x+2,n_out_start:n_out_end) = v_0(n_x+1,n_out_start:n_out_end);  % zero-gradient    
    
    % Advection-diffusion equation (predictor)
    [u_t, v_t] = Solving_Advection_Diffusion_2D(u_t, v_t, u_0, v_0, n_x, n_y, h_x, h_y, dt, nu);
    
    % Update boundary conditions for temporary velocity
    u_t(1,n_in_start:n_in_end)       = u_0(1,n_in_start:n_in_end);        % fixed velocity [m/s]
    u_t(n_x+1,n_out_start:n_out_end) = u_0(n_x+1,n_out_start:n_out_end);  % zero-gradient
    v_t(n_x+2,n_out_start:n_out_end) = v_0(n_x+2,n_out_start:n_out_end);  % zero-gradient 
    
    % Pressure equation
    [p_0, iter] = Poisson2D(p_0, u_t, v_t, alpha, n_x, n_y, h_x, h_y, dt, omega, max_iter, maximum_error, solver_option);
    
    % Correction on the velocity
    u_0(2:n_x,2:n_y+1) = u_t(2:n_x,2:n_y+1) - (dt/h_x)*(p_0(3:n_x+1,2:n_y+1)-p_0(2:n_x,2:n_y+1));
    v_0(2:n_x+1,2:n_y) = v_t(2:n_x+1,2:n_y) - (dt/h_y)*(p_0(2:n_x+1,3:n_y+1)-p_0(2:n_x+1,2:n_y));
    
    % Correction on outlet to ensure conservation of mass
    u_0(n_x+1,n_out_start:n_out_end)=u_t(n_x+1,n_out_start:n_out_end) - (dt/h_x)*(p_0(n_x+2,n_out_start:n_out_end)-p_0(n_x+1,n_out_start:n_out_end));
    
    %Numerical errors in equation solving may compromise mass conservation. Correcting outlet velocity is essential to ensure mass balance.
    Q_in  = mean(u_0(1,n_in_start:n_in_end))*A_in;         % inlet flow rate [m2/s]
    Q_out = mean(u_0(n_x+1,n_out_start:n_out_end))*A_out;  % outlet flow rate [m2/s]
    
    Q_in_values(is)   = Q_in;
    Q_out_values(is)  = Q_out;
    Q_diff            = Q_out - Q_in;
    Q_diff_values(is) = abs(Q_diff);
    
    % Update the plot
    semilogy(1:is, Q_diff_values(1:is)*10 , 'b.-', 'LineWidth', 1.3, 'markersize', 15);
    title('Mass Flow Rate Error vs Iteration');
    xlabel('Iteration');
    ylabel('% Mass Flow Rate Error [m^2/s]');
    set(gca, 'Box', 'on');
    set(gca, 'LineWidth', 1.3);
    box; drawnow;
    
    if (abs(Q_out)>1.e-6)
        u_0(n_x+1,n_out_start:n_out_end) = u_0(n_x+1,n_out_start:n_out_end)*abs(Q_in/Q_out);
    end
    
    % Print on the screen
    if (mod(is,50)==1)
        fprintf( 'Step: %d - Real Time: %f - iterations: %d - Mass Error: %f%%\n', is, t, iter, (Q_out-Q_in)/Q_in*10);
    end
    
    % Advance in time
    t = t+dt;
 
end

% Results nad post - processing                                                   %

% Field reconstruction
u_u(1:n_x+1,1:n_y+1) = 0.50*(u_0(1:n_x+1,2:n_y+2)+u_0(1:n_x+1,1:n_y+1));
v_v(1:n_x+1,1:n_y+1) = 0.50*(v_0(2:n_x+2,1:n_y+1)+v_0(1:n_x+1,1:n_y+1));
p_p(1:n_x+1,1:n_y+1) = 0.25*(p_0(1:n_x+1,1:n_y+1)+p_0(1:n_x+1,2:n_y+2) + p_0(2:n_x+2,1:n_y+1)+p_0(2:n_x+2,2:n_y+2));

line_colors = lines(3);

% Surface map: u-velocity
figure(1);
surf(X, Y, u_u', 'EdgeColor', 'none');
title('u-velocity'); 
xlabel('x'); 
ylabel('y');
shading interp; colorbar; colormap jet
axis([0 L_x 0 L_y]);
box; view(2)
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);

% Surface map: v-velocity
figure(2);
surf(X, Y, v_v', 'EdgeColor', 'none');
title('v-velocity');
xlabel('x'); 
ylabel('y');
shading interp; colorbar; colormap jet
axis([0 L_x 0 L_y]);
box; view(2)
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);

% Surface map: v-velocity
figure(3);
surf(X, Y, sqrt((u_u').^2+(v_v').^2), 'EdgeColor', 'none');
title('velocity magnitude');
xlabel('x'); 
ylabel('y');
shading interp; colorbar; colormap jet
axis([0 L_x 0 L_y]);
box; view(2)
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);

% Surface map: pressure
figure(4);
surf(X, Y, p_p', 'EdgeColor', 'none');
title('Pressure'); 
xlabel('x'); 
ylabel('y');
shading interp; colorbar; colormap jet
axis([0 L_x 0 L_y]);
box; view(2)
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);

% streamlines
figure(5);
sy = linspace(0, L_y, 100);
sx = L_x * 0.01 * ones(size(sy));
h = streamline(X, Y, u_u', v_v', sx, sy, [0.001, 100000]);
set(h, 'Color', 'k', 'LineWidth', 1.5);
hold on;
plot(sx, sy, 'ro', 'MarkerFaceColor', 'r');
axis([0 L_x 0 L_y]);
title('Streamlines');
xlabel('x');
ylabel('y');
hold off;
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);

% Surface map: velocity vectors
figure(6);
quiver(X, Y, u_u', v_v', 'color', 'b', 'LineWidth', 1.2, 'MaxHeadSize', 0.05, 'AutoScaleFactor', 2);
axis([0 L_x 0 L_y]);
title('Velocity Vector Field');
xlabel('x');
ylabel('y');
box;
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);

% Plot components along the vertical middle axis
figure(7);
plot(x, u_u(:, round(n_y/2)+1), 'LineWidth', 2.5, 'Color', line_colors(1, :));
hold on;
plot(x, v_v(:, round(n_y/2)+1), 'LineWidth', 2.5, 'Color', line_colors(2, :));
plot(x, sqrt((u_u(:, round(n_y/2)+1)).^2+(v_v(:, round(n_y/2)+1)).^2), 'LineWidth', 2.5, 'Color', line_colors(3, :));
hold off;
xlabel('x [m]');
ylabel('Velocity [m/s]');
title('Velocity along x-axis in y/L=0.5');
legend('x-velocity', 'y-velocity','velocity', 'Location', 'best');
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);
axis('square');

% Plot components along the horizontal middle axis
figure(8);
plot(u_u(round(n_x/2)+1,:),y, 'LineWidth', 2.5, 'Color', line_colors(1, :));
hold on;
plot(v_v(round(n_x/2)+1,:),y, 'LineWidth', 2.5, 'Color', line_colors(2, :));
plot(sqrt((u_u(round(n_x/2)+1,:)).^2+(v_v(round(n_x/2)+1,:)).^2), y, 'LineWidth', 2.5, 'Color', line_colors(3, :));
hold off;
xlabel('velocity [m/s]');
ylabel('y [m]');
title('velocity along y-axis in x/L=0.5');
legend('x-velocity', 'y-velocity', 'velocity', 'Location', 'best');
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);
axis('square');

% Plot pressure along the vertical and horizontal middle axis
figure(9);
plot(x, p_p(:, round(n_y/2)+1), 'LineWidth', 2.5, 'Color', line_colors(1, :));
xlabel('x [m]');
ylabel('Presure');
title('pressure along middle x-axis in L/2');
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);
axis('square');

% Plot pressure along the vertical and horizontal middle axis
figure(10);
plot(p_p(round(n_x/2)+1,:),y,   'LineWidth', 2.5, 'Color', line_colors(2, :));
xlabel('Presure');
ylabel('y [m]');
title('pressure along middle y-axis in L/2');
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.3);
axis('square');

fprintf( '\nPressure: %f | u-velocity: %f | v-velocity: %f \n\n', p_p(round(n_x/2)+1,round(n_y/2)+1),  u_u(round(n_x/2)+1,round(n_y/2)+1), v_v(round(n_x/2)+1,round(n_y/2)+1));         


% Poisson equation solver
function [p, iteration] = Poisson2D( p, u_t, v_t, alpha, n_x, n_y, h_x, h_y, dt, beta, max_iterations, max_error, solver_type)

    % SOR solver
    if (strcmp(solver_type,'SOR'))

        % Main loop
        for iteration=1:max_iterations
            
            for i=2:n_x+1
                for j=2:n_y+1

                    delta  = h_y/h_x*(p(i+1,j)+p(i-1,j))+h_x/h_y*(p(i,j+1)+p(i,j-1));
                    S      = 1/dt*( h_y*(u_t(i,j)-u_t(i-1,j)) + h_x*(v_t(i,j)-v_t(i,j-1)));
                    p(i,j) = beta*alpha(i,j)*( delta-S )+(1-beta)*p(i,j);

                end
            end
            
            % Estimate the error`
            epsilon     = 0.0;
            count_cells = 0;
            for i=2:n_x+1
                for j=2:n_y+1
    
                    delta   = h_y/h_x*(p(i+1,j)+p(i-1,j))+h_x/h_y*(p(i,j+1)+p(i,j-1));
                    S       = 1/dt*( h_y*(u_t(i,j)-u_t(i-1,j)) + h_x*(v_t(i,j)-v_t(i,j-1)));           
                    epsilon = epsilon+abs( p(i,j) - alpha(i,j)*( delta-S ) );
                    count_cells = count_cells + 1;

                end
            end
            epsilon = epsilon / count_cells;
            
            % Check the error
            if (epsilon <= max_error) % stop if converged
                break;
            end 
            
        end

    else    % Direct or GMRES solvers

        ne = (n_x+2)*(n_y+2);
        b = zeros(ne,1);

        % Fill main diagonal
        counter = 1;
        for i=1:ne
            I(counter) = i; J(counter) = i; V(counter) = 1.; counter = counter+1;
        end
    
        % Fill equations
        for i=2:n_x+1
            for j=2:n_y+1
                
                k = (n_x+2)*(j-1) + i;

                I(counter) = k; J(counter) = k+1;       V(counter) = -alpha(i,j)*h_y/h_x; counter = counter+1;
                I(counter) = k; J(counter) = k-1;       V(counter) = -alpha(i,j)*h_y/h_x; counter = counter+1;
                I(counter) = k; J(counter) = k+(n_x+2); V(counter) = -alpha(i,j)*h_x/h_y; counter = counter+1;
                I(counter) = k; J(counter) = k-(n_x+2); V(counter) = -alpha(i,j)*h_x/h_y; counter = counter+1;

                b(k) = -alpha(i,j)*(1/dt)*(h_y*(u_t(i,j)-u_t(i-1,j))+h_x*(v_t(i,j)-v_t(i,j-1)));
            
            end
        end
    
        M = sparse(I,J,V,ne,ne);
    
        if (strcmp(solver_type,'Direct'))
    
            p = M\b;
            iteration = 0;
    
        elseif (strcmp(solver_type,'GMRES'))
    
            tol = 1e-6;
            maxit = 10;
            [p,~,~,iteration] = gmres(M,b,[],tol,maxit,[],[],p(:));
            iteration=iteration(2);
    
        end
    
        p = reshape(p,[n_x+2 n_y+2]);
    
    end

end

% Advection-diffusion equation
function [u_t, v_t] = Solving_Advection_Diffusion_2D(u_t, v_t, u, v, n_x, n_y, h_x, h_y, dt, nu)
                       
    % solving u-velocity
    for i=2:n_x
        for j=2:n_y+1 
            
            u_e = (u(i,j)+u(i+1,j))/2;
            u_w = (u(i,j)+u(i-1,j))/2;
            u_n = u(i,j)+(u(i,j+1)-u(i,j))*0.50;
            u_s = u(i,j-1)+(u(i,j)-u(i,j-1))*0.50;
            v_n = v(i,j)+(v(i+1,j)-v(i,j))*0.50;
            v_s = v(i,j-1)+(v(i+1,j-1)-v(i,j-1))*0.50;
            
            u_e_2 = u_e^2 * h_y;
            u_w_2 = u_w^2 * h_y;
            u_n_v = u_n*v_n * h_x;
            u_s_v = u_s*v_s * h_x;
            
            Volume = h_x * h_y;
            Aria = (u_e_2-u_w_2+u_n_v-u_s_v)/Volume;
            
            D_e = nu*(u(i+1,j)-u(i,j))/h_x*h_y;
            D_w = nu*(u(i,j)-u(i-1,j))/h_x*h_y;
            D_n = nu*(u(i,j+1)-u(i,j))/h_y*h_x;
            D_s = nu*(u(i,j)-u(i,j-1))/h_y*h_x;
            D_final = (D_e-D_w+D_n-D_s)/Volume;
            
            u_t(i,j)=u(i,j)+dt*(-Aria+D_final);
            
        end
    end
    
    % solving v-velocity
    for i=2:n_x+1
        for j=2:n_y 
            
            v_n = (v(i,j)+v(i,j+1))/2;
            v_s = (v(i,j)+v(i,j-1))/2;
            v_e = v(i,j)+(v(i+1,j)-v(i,j))*0.50;
            v_w = v(i-1,j)+(v(i,j)-v(i-1,j))*0.50;
            u_e = u(i,j)+(u(i,j+1)-u(i,j))*0.50;
            u_w = u(i-1,j)+(u(i-1,j+1)-u(i-1,j))*0.50;
            
            v_n_2 = v_n^2 * h_x;
            v_s_2 = v_s^2 * h_x;
            v_e_u = v_e*u_e * h_y;
            v_w_u = v_w*u_w * h_y;
            
            Volume = h_x * h_y;
            Aria = (v_n_2 - v_s_2 + v_e_u - v_w_u)/Volume;
            
            D_e = nu*(v(i+1,j)-v(i,j))/h_x*h_y;
            D_w = nu*(v(i,j)-v(i-1,j))/h_x*h_y;
            D_n = nu*(v(i,j+1)-v(i,j))/h_y*h_x;
            D_s = nu*(v(i,j)-v(i,j-1))/h_y*h_x;
            D_final = (D_e-D_w+D_n-D_s)/Volume;
            
            v_t(i,j)=v(i,j)+dt*(-Aria+D_final);
            
        end
    end
    
end
