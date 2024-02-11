% Solving Compression Wave Flow Field Numerically %
% Ali Ramezani 401126144%
% Advance CFD %

clc
clear
close all

% Initial Conditions: ISA Sea Level Standards
M_infinity     = 2.93;
% 2.5-2.93
T_infinity     = 298;
p_infinity     = 101325;
R_air_1        = 287;
gamma_heat_air = 1.4;
rho_infinity   = p_infinity/(R_air_1*T_infinity);
a_infinity     = sqrt(gamma_heat_air*R_air_1*T_infinity);
V_infinity     = M_infinity*a_infinity;

% Definition of Physical Domain %
H_physical  = 1.0;
L_physical  = 1.0;
E_compress  = 0.8;
angle_12    = 35;
Theta_angel = deg2rad(27.416);
% 18.829 - 27.416
x_initial   = 0;

% Parameters for Grid Size %
J_n                 = 100;
delta_eta_incerment = 1/(J_n - 1);
iteration_index      = 1;
maximum_iterations  = 1e5;
CFL_number          = 0.5;
alpha_visco         = 1;

% Definition of Computational/Transformed Domain %
eta_plane = (0:delta_eta_incerment:1)';

% Initial data line includes:
y_s_wall         = 0;
h_1              = H_physical;
y                = eta_plane*h_1 + y_s_wall;
deta_dx(1:J_n,1) = 0;
deta_dy          = 1/h_1;

% Initialization of Primitive Flow Variables %
% Defining Primitive Flow Variables at Initial Data Line based on Free-Stream Conditions %
p(1:J_n,1)   = p_infinity;
rho(1:J_n,1) = rho_infinity;
u(1:J_n,1)   = V_infinity;
v(1:J_n,1)   = 0;

% Initializing Flux Terms from Primitive Flow Variables %
F_1 = rho.*u;
F_2 = rho.*u.^2 + p;
F_3 = rho.*u.*v;
F_4 = gamma_heat_air/(gamma_heat_air - 1)*p.*u + rho.*u.*(u.^2 + v.^2)/2;
G_1 = rho.*v;
G_2 = F_3;
G_3 = rho.*v.^2 + p;
G_4 = gamma_heat_air/(gamma_heat_air - 1)*p.*v + rho.*v.*(u.^2 + v.^2)/2;

% Implementing Finite Difference Solver %
% Using MacCormack's Predictor-Corrector Scheme for a straightforward solution.
SF1   = zeros(J_n,1);
SF2   = zeros(J_n,1);
SF3   = zeros(J_n,1);
SF4   = zeros(J_n,1);
SF1_p = zeros(J_n,1);
SF2_p = zeros(J_n,1);
SF3_p = zeros(J_n,1);
SF4_p = zeros(J_n,1);

while (x_initial(iteration_index) < L_physical && iteration_index < maximum_iterations)

    M_a(:,iteration_index)   = sqrt(u(:,iteration_index).^2 + v(:,iteration_index).^2)./sqrt(gamma_heat_air*p(:,iteration_index)./rho(:,iteration_index));
    le_running              = abs(tan(atan(v(:,iteration_index)./u(:,iteration_index)) + asin(1./M_a(:,iteration_index))));
    ri_running              = abs(tan(atan(v(:,iteration_index)./u(:,iteration_index)) - asin(1./M_a(:,iteration_index))));
    delta_x(iteration_index) = CFL_number*delta_eta_incerment*h_1(iteration_index)/max(max(le_running),max(ri_running));
    
    for j = 1:J_n-1
        dF1_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F_1(j+1,iteration_index) - F_1(j,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G_1(j+1,iteration_index) - G_1(j,iteration_index))/delta_eta_incerment;
        dF2_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F_2(j+1,iteration_index) - F_2(j,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G_2(j+1,iteration_index) - G_2(j,iteration_index))/delta_eta_incerment;
        dF3_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F_3(j+1,iteration_index) - F_3(j,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G_3(j+1,iteration_index) - G_3(j,iteration_index))/delta_eta_incerment;
        dF4_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F_4(j+1,iteration_index) - F_4(j,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G_4(j+1,iteration_index) - G_4(j,iteration_index))/delta_eta_incerment;
    end
    
    j = J_n;
    dF1_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F_1(j,iteration_index) - F_1(j-1,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G_1(j,iteration_index) - G_1(j-1,iteration_index))/delta_eta_incerment;
    dF2_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F_2(j,iteration_index) - F_2(j-1,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G_2(j,iteration_index) - G_2(j-1,iteration_index))/delta_eta_incerment;
    dF3_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F_3(j,iteration_index) - F_3(j-1,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G_3(j,iteration_index) - G_3(j-1,iteration_index))/delta_eta_incerment;
    dF4_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F_4(j,iteration_index) - F_4(j-1,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G_4(j,iteration_index) - G_4(j-1,iteration_index))/delta_eta_incerment;
    
    for j = 2:J_n-1
        p_sensor = alpha_visco*abs(p(j+1,iteration_index) - 2*p(j,iteration_index) + p(j-1,iteration_index))/(p(j+1,iteration_index) + 2*p(j,iteration_index) + p(j-1,iteration_index));
        SF1(j) = p_sensor*(F_1(j+1,iteration_index) - 2*F_1(j,iteration_index) + F_1(j-1,iteration_index));
        SF2(j) = p_sensor*(F_2(j+1,iteration_index) - 2*F_2(j,iteration_index) + F_2(j-1,iteration_index));
        SF3(j) = p_sensor*(F_3(j+1,iteration_index) - 2*F_3(j,iteration_index) + F_3(j-1,iteration_index));
        SF4(j) = p_sensor*(F_4(j+1,iteration_index) - 2*F_4(j,iteration_index) + F_4(j-1,iteration_index));
    end
    
    F1_p(:,iteration_index) = F_1(:,iteration_index) + dF1_dx(:,iteration_index)*delta_x(iteration_index) + SF1;
    F2_p(:,iteration_index) = F_2(:,iteration_index) + dF2_dx(:,iteration_index)*delta_x(iteration_index) + SF2;
    F3_p(:,iteration_index) = F_3(:,iteration_index) + dF3_dx(:,iteration_index)*delta_x(iteration_index) + SF3;
    F4_p(:,iteration_index) = F_4(:,iteration_index) + dF4_dx(:,iteration_index)*delta_x(iteration_index) + SF4;
    
    [G1_p(:,iteration_index), G2_p(:,iteration_index), G3_p(:,iteration_index), G4_p(:,iteration_index)] = F2G(F1_p(:,iteration_index), F2_p(:,iteration_index), F3_p(:,iteration_index), F4_p(:,iteration_index), gamma_heat_air);
    p_p(:,iteration_index) = F2Primitive(F1_p(:,iteration_index), F2_p(:,iteration_index), F3_p(:,iteration_index), F4_p(:,iteration_index), gamma_heat_air);

    for j = 2:J_n
        dF1_p_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F1_p(j,iteration_index) - F1_p(j-1,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G1_p(j,iteration_index) - G1_p(j-1,iteration_index))/delta_eta_incerment;
        dF2_p_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F2_p(j,iteration_index) - F2_p(j-1,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G2_p(j,iteration_index) - G2_p(j-1,iteration_index))/delta_eta_incerment;
        dF3_p_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F3_p(j,iteration_index) - F3_p(j-1,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G3_p(j,iteration_index) - G3_p(j-1,iteration_index))/delta_eta_incerment;
        dF4_p_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F4_p(j,iteration_index) - F4_p(j-1,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G4_p(j,iteration_index) - G4_p(j-1,iteration_index))/delta_eta_incerment;
    end
    
    j = 1;
    dF1_p_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F1_p(j+1,iteration_index) - F1_p(j,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G1_p(j+1,iteration_index) - G1_p(j,iteration_index))/delta_eta_incerment;
    dF2_p_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F2_p(j+1,iteration_index) - F2_p(j,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G2_p(j+1,iteration_index) - G2_p(j,iteration_index))/delta_eta_incerment;
    dF3_p_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F3_p(j+1,iteration_index) - F3_p(j,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G3_p(j+1,iteration_index) - G3_p(j,iteration_index))/delta_eta_incerment;
    dF4_p_dx(j,iteration_index) = - deta_dx(j,iteration_index)*(F4_p(j+1,iteration_index) - F4_p(j,iteration_index))/delta_eta_incerment - deta_dy(iteration_index)*(G4_p(j+1,iteration_index) - G4_p(j,iteration_index))/delta_eta_incerment;
    
    for j = 2:J_n-1
        p_p_sensor = alpha_visco*abs(p_p(j+1,iteration_index) - 2*p_p(j,iteration_index) + p_p(j-1,iteration_index))/(p_p(j+1,iteration_index) + 2*p_p(j,iteration_index) + p_p(j-1,iteration_index));
        SF1_p(j) = p_p_sensor*(F1_p(j+1,iteration_index) - 2*F1_p(j,iteration_index) + F1_p(j-1,iteration_index));
        SF2_p(j) = p_p_sensor*(F2_p(j+1,iteration_index) - 2*F2_p(j,iteration_index) + F2_p(j-1,iteration_index));
        SF3_p(j) = p_p_sensor*(F3_p(j+1,iteration_index) - 2*F3_p(j,iteration_index) + F3_p(j-1,iteration_index));
        SF4_p(j) = p_p_sensor*(F4_p(j+1,iteration_index) - 2*F4_p(j,iteration_index) + F4_p(j-1,iteration_index));
    end
    
    F_1(:,iteration_index+1) = F_1(:,iteration_index) + 1/2*(dF1_dx(:,iteration_index) + dF1_p_dx(:,iteration_index))*delta_x(iteration_index) + SF1_p;
    F_2(:,iteration_index+1) = F_2(:,iteration_index) + 1/2*(dF2_dx(:,iteration_index) + dF2_p_dx(:,iteration_index))*delta_x(iteration_index) + SF2_p;
    F_3(:,iteration_index+1) = F_3(:,iteration_index) + 1/2*(dF3_dx(:,iteration_index) + dF3_p_dx(:,iteration_index))*delta_x(iteration_index) + SF3_p;
    F_4(:,iteration_index+1) = F_4(:,iteration_index) + 1/2*(dF4_dx(:,iteration_index) + dF4_p_dx(:,iteration_index))*delta_x(iteration_index) + SF4_p;
    
    [p(2:J_n,iteration_index+1), rho(2:J_n,iteration_index+1), u(2:J_n,iteration_index+1), v(2:J_n,iteration_index+1)] = F2Primitive(F_1(2:J_n,iteration_index+1), F_2(2:J_n,iteration_index+1), F_3(2:J_n,iteration_index+1), F_4(2:J_n,iteration_index+1), gamma_heat_air);
    
    [p_old, rho_old, u_old, v_old] = F2Primitive(F_1(1,iteration_index+1), F_2(1,iteration_index+1), F_3(1,iteration_index+1), F_4(1,iteration_index+1), gamma_heat_air);
    M_old = sqrt(u_old^2 + v_old^2)/sqrt(gamma_heat_air*p_old/rho_old);
    
    x_initial(iteration_index+1) = x_initial(iteration_index) + delta_x(iteration_index);
    if (x_initial(iteration_index+1) <= E_compress)
        y_s_wall(iteration_index+1) = 0;
        h_1(iteration_index+1) = H_physical;
        theta = 0;
        deta_dx(:,iteration_index+1) = 0;
    else
        y_s_wall(iteration_index+1) = (x_initial(iteration_index+1) - x_initial(find(x_initial >= E_compress,1)))*tan(Theta_angel); % Since we can't control if the there will be a grid point exactly at x = E.
        h_1(iteration_index+1) = H_physical - (x_initial(iteration_index+1) - x_initial(find(x_initial >= E_compress,1)))*tan(Theta_angel);
        theta =  Theta_angel;
        deta_dx(:,iteration_index+1) = (eta_plane - 1)*tan(Theta_angel)/h_1(iteration_index+1);
    end
    y(:,iteration_index+1) = eta_plane*h_1(iteration_index+1) + y_s_wall(iteration_index+1);
    
    phi = atan(v_old/u_old) - theta;
    nu_old = sqrt((gamma_heat_air + 1)/(gamma_heat_air - 1))*atan(sqrt((gamma_heat_air - 1)*(M_old^2 - 1)/(gamma_heat_air + 1))) - atan(sqrt(M_old^2 - 1));
    nu_new = nu_old + phi;
    
    M_new = find_PM_Mach(gamma_heat_air,nu_new,[1.01,M_old + 1],1e-8);
    
    p(1,iteration_index+1)   = p_old*((1 + (gamma_heat_air - 1)/2*M_old^2)/(1 + (gamma_heat_air - 1)/2*M_new^2))^(gamma_heat_air/(gamma_heat_air - 1));
    rho(1,iteration_index+1) = rho_old*((1 + (gamma_heat_air - 1)/2*M_old^2)/(1 + (gamma_heat_air - 1)/2*M_new^2))^(1/(gamma_heat_air - 1));
    V_t                      = u_old.*cos(theta) + v_old.*sin(theta);
    u(1,iteration_index+1)   = V_t.*cos(theta);
    v(1,iteration_index+1)   = V_t.*sin(theta);
    
    F_1(1,iteration_index+1) = rho(1,iteration_index+1).*u(1,iteration_index+1);
    F_2(1,iteration_index+1) = rho(1,iteration_index+1).*u(1,iteration_index+1).^2 + p(1,iteration_index+1);
    F_3(1,iteration_index+1) = rho(1,iteration_index+1).*u(1,iteration_index+1).*v(1,iteration_index+1);
    F_4(1,iteration_index+1) = gamma_heat_air/(gamma_heat_air - 1)*p(1,iteration_index+1).*u(1,iteration_index+1) + rho(1,iteration_index+1).*u(1,iteration_index+1).*(u(1,iteration_index+1).^2 + v(1,iteration_index+1).^2)/2;
    
    deta_dy(iteration_index+1) = 1/h_1(iteration_index+1);
    
    [G_1(:,iteration_index+1), G_2(:,iteration_index+1), G_3(:,iteration_index+1), G_4(:,iteration_index+1)] = F2G(F_1(:,iteration_index+1), F_2(:,iteration_index+1), F_3(:,iteration_index+1), F_4(:,iteration_index+1), gamma_heat_air);
    
    iteration_index = iteration_index + 1
    
    %error_1 = G1(iteration)-G1(iteration-1);
    %error_2 = G1(iteration)-G1(iteration-1);
    %error_3 = G1(iteration)-G1(iteration-1);
    %error_4 = G1(iteration)-G1(iteration-1);
    %semilogy(iteration,error_1) 
    %hold on
    %semilogy(iteration,error_2)
    %hold on
    %semilogy(iteration,error_3)
    %hold on    
    %semilogy(iteration,error_4) 
    %hold on
    %xlabel('Iteration')
    %ylabel('Residual')
end

if (iteration_index == maximum_iterations)
    error('The maximum number of iterations has been reached. The simulation is not completed')
end

% Other fundamental flow-field properties
Tem  = p./(R_air_1*rho);
Valu = sqrt(u.^2 + v.^2);
M_a  = Valu./sqrt(gamma_heat_air*R_air_1*Tem);

% Counntors and Plots
degreeSymbol_a = char(176);
M_infinity     = 2.5;
% Pressure
figure(1)
surf(x_initial,y,(p*10^(-3)))
shading interp; colorbar; colormap jet
title(['Pressure [kPa] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x')
ylabel('y')
axis image; box; view(2)

% Density
figure(2)
surf(x_initial,y,rho)
shading interp; colorbar; colormap jet
title(['Density [kg/m^3] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x')
ylabel('y')
axis image; box; view(2)

% x velocity
figure(3)
surf(x_initial,y,u)
shading interp; colorbar; colormap jet
title(['Velocity x [m/s] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x')
ylabel('y')
axis image; box; view(2)

% y velocity
figure(4)
surf(x_initial,y,v)
shading interp; colorbar; colormap jet
title(['Velocity y [m/s] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x')
ylabel('y')
axis image; box; view(2)

% velocity magnitude
figure(5)
surf(x_initial,y,Valu)
shading interp; colorbar; colormap jet
title(['Velocity magnitude [m/s] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x')
ylabel('y')
axis image; box; view(2)

% Temperature
figure(6)
surf(x_initial,y,Tem)
shading interp; colorbar; colormap jet
title(['Temperature [K] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x')
ylabel('y')
axis image; box; view(2)

% Mach number
figure(7)
surf(x_initial,y,M_a)
shading interp; colorbar; colormap jet
title(['Mach number (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x')
ylabel('y')
axis image; box; view(2)

% Non-Dimensional Countors

% Nondimensional Pressure
figure(8)
surf(x_initial/L_physical, y/L_physical, p/p_infinity)
shading interp; colorbar; colormap jet
title(['Non-dimensional pressure (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x/L')
ylabel('y/L')
axis image; box; view(2)

% Nondimensional Density
figure(9)
surf(x_initial/L_physical, y/L_physical, rho/rho_infinity)
shading interp; colorbar; colormap jet
title(['Non-dimensional Density [kg/m^3] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x/L')
ylabel('y/L')
axis image; box; view(2)

% Nondimensional x Velocity
figure(10)
surf(x_initial/L_physical, y/L_physical, u/V_infinity)
shading interp; colorbar; colormap jet
title(['Non-dimensional x Velocity [m/s] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x/L')
ylabel('y/L')
axis image; box; view(2)

% Nondimensional y Velocity
figure(11)
surf(x_initial/L_physical, y/L_physical, v/max(v(:)))
shading interp; colorbar; colormap jet
title(['Non-dimensional y Velocity [m/s] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x/L')
ylabel('y/L')
axis image; box; view(2)

% Nondimensional Velocity
figure(12)
surf(x_initial/L_physical, y/L_physical, v/max(Valu(:)))
shading interp; colorbar; colormap jet
title(['Non-dimensional Velocity [m/s] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x/L')
ylabel('y/L')
axis image; box; view(2)

% Nondimensional Temperature
figure(13)
surf(x_initial/L_physical, y/L_physical, Tem/T_infinity)
shading interp; colorbar; colormap jet
title(['Non-dimensional Temperature [K] (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x/L')
ylabel('y/L')
axis image; box; view(2)

% Nondimensional Mach 
figure(14)
surf(x_initial/L_physical, y/L_physical, M_a/M_infinity)
shading interp; colorbar; colormap jet
title(['Non-dimensional Mach number (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
xlabel('x/L')
ylabel('y/L')
axis image; box; view(2)

% Pressure Plot in x/L = 1
figure(15)
plot((p(:,end))*(10^(-3)),y(:, end)/L_physical,'Linewidth',3, 'Color', '[0.5 0.2 0.8]')
grid on
xlabel('Pressure [kPa] (in x/L = 1)')
ylabel('y/L')
title(['Freestream conditions: (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
ylim([0 1])

figure(16)
plot(rho(:,end),y(:, end)/L_physical,'Linewidth',3, 'Color', '[0.5 0.2 0.8]')
grid on
xlabel('Density [kg/m^3] (in x/L = 1)')
ylabel('y/L')
title(['Freestream conditions: (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
ylim([0 1])

figure(17)
plot(u(:,end),y(:, end)/L_physical,'Linewidth',3, 'Color', '[0.5 0.2 0.8]')
grid on
xlabel('x Velocity [m/s] (in x/L = 1)')
ylabel('y/L')
title(['Freestream conditions: (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
ylim([0 1])

figure(18)
plot(abs(v(:,end)),y(:, end)/L_physical,'Linewidth',3, 'Color', '[0.5 0.2 0.8]')
grid on
xlabel('y Velocity [m/s] (in x/L = 1)')
ylabel('y/L')
title(['Freestream conditions: (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
ylim([0 1])

figure(19)
plot(Valu(:,end),y(:, end)/L_physical,'Linewidth',3, 'Color', '[0.5 0.2 0.8]')
grid on
xlabel('Velocity Magnitude [m/s]  (in x/L = 1)')
ylabel('y/L')
title(['Freestream conditions: (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
ylim([0 1])

figure(20)
plot(Tem(:,end),y(:, end)/L_physical,'Linewidth',3, 'Color', '[0.5 0.2 0.8]')
grid on
xlabel('Temperature [K]  (in x/L = 1)')
ylabel('y/L')
title(['Freestream conditions: (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
ylim([0 1])

figure(21)
plot(M_a(:,end),y(:, end)/L_physical,'Linewidth',3, 'Color', '[0.5 0.2 0.8]')
grid on
xlabel('Mach Number (in x/L = 1)')
ylabel('y/L')
title(['Freestream conditions: (M_{\infty} = ', num2str(M_infinity), ' | CFL = ', num2str(CFL_number),' | \gamma = ',num2str(gamma_heat_air), ' | \theta = ', num2str(angle_12), degreeSymbol_a, ')'])
ylim([0 1])


% Calculate maximum and minimum values
max_pressure = max(p(:));     min_p   = min(p(:));
max_rho      = max(rho(:));   min_rho = min(rho(:));
max_u        = max(u(:));     min_u   = min(u(:));
max_val      = max(v(:));     min_v   = min(v(:));
max_vm       = max(Valu(:));  min_vm  = min(Valu(:));
max_Tem      = max(Tem(:));   min_T   = min(Tem(:));
max_Ma       = max(M_a(:));   min_M   = min(M_a(:));

% Create a table
data = table({max_pressure; min_p; max_rho; min_rho; max_u; min_u; max_val; min_v; max_vm; min_vm; max_Tem; min_T; max_Ma; min_M}, ...
             'VariableNames', {'Values'},'RowNames', {'Max Pressure', 'Min Pressure', 'Max Density', 'Min Density', 'Max x Velocity', 'Min x Velocity', 'Max y Velocity', 'Min y Velocity', 'Max Velocity Magnitude', 'Min Velocity Magnitude', 'Max Temperature', 'Min Temperature', 'Max Mach Number', 'Min Mach Number'});

% Display the table
disp(data);

%Derive primitive variables from flux vector F in a 2D, inviscid, calorically perfect gas using a specific function.
function [p, rho, u, v] = F2Primitive(F1, F2, F3, F4, gamma)
    A = F3.^2./(2*F1) - F4;
    B = gamma/(gamma - 1)*F1.*F2;
    C = - (gamma + 1)/(2*(gamma - 1))*F1.^3;
    rho = (- B + sqrt(B.^2 - 4*A.*C))./(2*A);
    u = F1./rho;
    v = F3./F1;
    p = F2 - F1.^2./rho;
end

%Compute flux vector G from solution vector F in a 2D, inviscid, calorically perfect gas using a specific function.
function [G_1, G_2, G_3, G_4] = F2G(F1, F2, F3, F4, gamma)
    A_1 = F3.^2./(2*F1) - F4;
    B_1 = gamma/(gamma - 1)*F1.*F2;
    C_1 = - (gamma + 1)/(2*(gamma - 1))*F1.^3;
    rho = (- B_1 + sqrt(B_1.^2 - 4*A_1.*C_1))./(2*A_1);

    G_1 = rho.*F3./F1;
    G_2 = F3;
    G_3 = rho.*(F3./F1).^2 + F2 - F1.^2./rho;
    G_4 = gamma/(gamma - 1)*(F2 - F1.^2./rho).*F3./F1 + rho/2.*F3./F1.*((F1./rho).^2 + (F3./F1).^2);
end

% Create the find_PM_Mach function
function M_1 = find_PM_Mach(gamma, nu_new, zero_int, tol)
    M_maximum  = 50;
    max_it     = 1e3;
    PM_fun     = @(M) sqrt((gamma + 1)/(gamma - 1))*atan(sqrt((gamma - 1)*(M^2 - 1)/(gamma + 1))) - atan(sqrt(M^2 - 1));
    nu_maximum = PM_fun(M_maximum);

    if gamma < 1
        error('The ratio of specific heats cannot be smaller than 1.')
    elseif gamma > 2
        error('The ratio of specific heats is out of range. The maximum allowed value is 2.')
    end
    if nu_new < 0
        error('Prandtl-Meyer function value cannot be negative.')
    elseif nu_new > nu_maximum
        error('Prandtl-Meyer function value is out of range. Maximum allowed value corresponds to a Mach number of 50.')
    elseif nu_new == 0
        M_1 = 1;
        return
    elseif nu_new == nu_maximum
        M_1 = 50;
        return
    end
    if min(zero_int) <= 1
        error('The solution interval has to contain Mach numbers larger than 1.')
    end
    if tol < eps
        error('The tolerance value has to be larger than the floating-point relative accuracy (given by eps).')
    end

    zero_fun = @(M) PM_fun(M) - nu_new;
    if zero_fun(zero_int(1))*zero_fun(zero_int(2)) > 0
        error('The function sign at the extremes of the solution interval has to be different.')
    end

    M_root = (zero_int(1) + zero_int(2))/2;
    i_count = 1;
    while abs(zero_fun(M_root)) > tol && i_count < max_it
        if zero_fun(zero_int(1))*zero_fun(M_root) < 0
            zero_int(2) = M_root;
        else
            zero_int(1) = M_root;
        end
        M_root = (zero_int(1) + zero_int(2))/2;
        i_count = i_count + 1;
    end

    if i_count == max_it
        error('Bisection algorithm did not find the root.')
    end

    M_1 = M_root;
end
