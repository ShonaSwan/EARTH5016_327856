%Session 5
%Two-dimensional (2D) model of heat transport by advection and diffusion.
%Explicit Time Integration, 


%***** RUN 2D ADVECTION DIFFUSION MODEL (Explicit) *******

% clear workspace
clear all; close all; %clc;

% set model parameters
W  = 500;          % domain width [m]
L =  500;          % domain length [m]
w  = 50;           % magma width [m]
Nx  = 200;          % grid size x direction
Nz  = 200;          % grid size z direction
N = 200
dx = W/Nx;          % grid spacing in x direction 
dz = L/Nz           % grid spacing in z direction   

T0   = 100;        % initial background temperature [C]
dT   = 1000;       % initial temperature peak amplitude [C]
wT   = 20;         % initial temperature peak width [m]
k0   = 1e-6;       % heat diffusivity [m2/s]
u0   = 1e-6;       % advection speed in  x-direction[m/s]
v0   = 1e-6;       % advection speed in z-direction [m/s]
BC   = 'periodic'; % boundary condition option flag ('insulating', 'periodic')
ADVN = 'WENO5';    % advection scheme ('UPW1', 'CFD2', 'UPW3', 'WENO5')
TINT = 'FE1';      % time integration scheme ('FE1','RK2','RK4')

yr    = 3600*24*365;  % seconds per year [s]
tend  = W/u0;         % stopping time [s]
CFL   = 1/16;         % Time step limiter
nop   = 100;          % output figure produced every 'nop' steps

%*****  RUN MODEL

%*****  Initialise Model Setup

% create x-coordinate vectors
xc = dx/2:dx:W-dx/2;    % coordinate vector for cell centre positions  x direction [m]
zc = dz/2:dz:L-dz/2;    % coordinate vector for cell centre positions  z direction [m]
xf = 0:dx:W;            % coordinate vector for cell face positions  x direction [m]
zf = 0:dz:L;            % coordinate vector for cell face positions  z direction [m]
[xc, zc] = meshgrid(xc, zc);

% set time step size
dt = CFL * min([dx/u0, dz/v0, dx^2/k0, dz^2/k0]); % time step [s]

% set up index array for boundary conditions
switch BC
    case 'periodic'
        % example periodic indexing for N=4
        %   [4,1,2,3,4,1]    % 3-point stencil
        % [3,4,1,2,3,4,1,2]  % 5-point stencil
        ind3_x = [        Nx,1:Nx,1    ];
        ind5_x = [    Nx-1,Nx,1:Nx,1,2  ];
        ind7_x = [Nx-2,Nx-1,Nx,1:Nx,1,2,3];

        % example periodic indexing for N=4
        %   [4,1,2,3,4,1]    % 3-point stencil
        % [3,4,1,2,3,4,1,2]  % 5-point stencil
        ind3_z = [        Nz,1:Nz,1    ];
        ind5_z = [    Nz-1,Nz,1:Nz,1,2  ];
        ind7_z = [Nz-2,Nz-1,Nz,1:Nz,1,2,3];

    case 'insulating'
        % example non-periodic indexing for N=4
        %   [1,1,2,3,4,4]      % 3-point stencil
        % [1,1,1,2,3,4,4,4]    % 5-point stencil
        ind3_x = [    1,1:Nx,Nx    ];
        ind5_x = [  1,1,1:Nx,Nx,Nx  ];
        ind7_x = [1,1,1,1:Nx,Nx,Nx,Nx];

        % example non-periodic indexing for N=4
        %   [1,1,2,3,4,4]      % 3-point stencil
        % [1,1,1,2,3,4,4,4]    % 5-point stencil
        ind3_z = [    1,1:Nz,Nz    ];
        ind5_z = [  1,1,1:Nz,Nz,Nz  ];
        ind7_z  = [1,1,1,1:Nz,Nz,Nz,Nz];

end

% set initial condition for temperature at cell centres
T   = T0 + dT .* exp(-((xc-W/2).^2 + (zc - L/2).^2)/(4*wT^2));     % initialise T array at Tr
Tin = T;                                         % store initial condition
Ta  = T;                                         % initialise analytical solution

% initialise output figure with initial condition

% figure(1); clf
% makefig2D(xc,zc,T,Tin,0)

figure(2); clf
makefig3D(xc,zc,T,Tin,0)

%*****  Solve Model Equations

t = 0;  % initial time [s]
k = 0;  % initial time step count

while t <= tend

    % increment time and step count
    t = t+dt;
    k = k+1;

    switch TINT
        case 'FE1'  % 1st-order Forward Euler time integration scheme
            
            dTdt = diffusion2D(T,k0,dx,dz,ind3_x,ind3_z) + advection2D(T,u0,v0,dx,dz,ind7_x,ind7_z,ADVN);

            T = T + dTdt * dt;

        case 'RK2'  % 2nd-order Runge-Kutta time integration scheme
            
            dTdt1 = diffusion2D(T           ,k0,dx,dz,ind3_x,ind3_z) + advection2D(T           ,u0,v0,dx,dz,ind7_x,ind7_z,ADVN);
            dTdt2 = diffusion2D(T+dTdt1*dt/2,k0,dx,dz,ind3_x,ind3_z) + advection2D(T+dTdt1*dt/2,u0,v0,dx,dz,ind7_x,ind7_z,ADVN);

            T = T + dTdt2 * dt;

        case 'HE2'  % 2nd-order Heun's time integration scheme
            
            dTdt1 = diffusion2D(T         ,k0,dx,dz,ind3_x,ind3_z) + advection2D(T         ,u0,v0,dx,dz,ind7_x,ind7_z,ADVN);
            dTdt2 = diffusion2D(T+dTdt1*dt,k0,dx,dz,ind3_x,ind3_z) + advection2D(T+dTdt1*dt,u0,v0,dx,dz,ind7_x,ind7_z,ADVN);

            T = T + (dTdt1 + dTdt2)/2 * dt;

        case 'RK4'  % 4th-order Runge-Kutta time integration scheme
            
            dTdt1 = diffusion2D(T           ,k0,dx,dz,ind3_x,ind3_z) + advection2D(T           ,u0,v0,dx,dz,ind7_x,ind7_z,ADVN);
            dTdt2 = diffusion2D(T+dTdt1/2*dt,k0,dx,dz,ind3_x,ind3_z) + advection2D(T+dTdt1/2*dt,u0,v0,dx,dz,ind7_x,ind7_z,ADVN);
            dTdt3 = diffusion2D(T+dTdt2/2*dt,k0,dx,dz,ind3_x,ind3_z) + advection2D(T+dTdt2/2*dt,u0,v0,dx,dz,ind7_x,ind7_z,ADVN);
            dTdt4 = diffusion2D(T+dTdt3  *dt,k0,dx,dz,ind3_x,ind3_z) + advection2D(T+dTdt3  *dt,u0,v0,dx,dz,ind7_x,ind7_z,ADVN);

            T = T + (dTdt1 + 2*dTdt2 + 2*dTdt3 + dTdt4)/6 * dt;

    end

    % get analytical solution at time t
    wTt = sqrt(wT^2 + 1*k0*t);
    Ta = T0 + dT * (wT^2 / wTt^2) * ( ...
    exp(-(((xc - W / 2 - u0 * t).^2) + ((zc - L / 2 - v0 * t).^2)) / (4 * wTt^2)) + ...
    exp(-(((xc + W / 2 - u0 * t).^2) + ((zc + L / 2 - v0 * t).^2)) / (4 * wTt^2)));

    % plot model progress
    if ~mod(k,nop)
        % figure(1)
        % makefig2D(xc,zc,T,Tin,t/yr);
        figure(2)
        makefig3D(xc,zc,T,Tin,t/yr);
    end

end

%*****  calculate numerical error norm
Err = norm(T(:) - Ta(:),2)./norm(Ta(:),2);
disp(' ');
disp(['Advection scheme: ',ADVN]);
disp(['Time integration scheme: ',TINT]);
disp(['Numerical error = ',num2str(Err)]);
disp(' ');

%*****  Utility Functions  ************************************************

% Function to make output figure
% function makefig2D(x,z,T,Tin,t)
% figure(1); 
% contourf(x, z, T, 50, 'LineStyle', 'none'); colorbar;
%     xlabel('x [m]', 'FontSize', 15);
%     ylabel('z [m]', 'FontSize', 15);
%     title(['Temperature; time = ', num2str(t), ' yr'], 'FontSize', 18);
%     drawnow;
% 
% end

% Function to make output figure in 3D
function makefig3D(x, z, T, Tin, t)
    figure(2);
    surf(x, z, T, 'EdgeColor', 'none'); % Create 3D surface plot
    colormap(jet); % Set color map
    colorbar; % Display color bar
    zlim([0 1200]);
    xlabel('x [m]', 'FontSize', 15);
    ylabel('z [m]', 'FontSize', 15);
    zlabel('Temperature [°C]', 'FontSize', 15);
    title(['Temperature; time = ', num2str(t), ' yr'], 'FontSize', 18);
    view(45, 45); % Adjust view angle for better visualization
    shading interp; % Smooth the surface
    drawnow;
end


% Function to calculate diffusion rate
function dTdt = diffusion2D(f,k0,dx,dz,ind_x,ind_z)
 % x-direction diffusion
    fx = f(:, ind_x); % Extend grid in the x-direction using stencil indices
    qx = -k0 .* diff(fx, 1, 2) / dx; % Heat flux in x-direction
    dqx = diff(qx, 1, 2) / dx; % Divergence of flux in x-direction

    % z-direction diffusion
    fz = f(ind_z, :); % Extend grid in the z-direction using stencil indices
    qz = -k0 .* diff(fz, 1, 1) / dz; % Heat flux in z-direction
    dqz = diff(qz, 1, 1) / dz; % Divergence of flux in z-direction

    % Combine contributions from both directions
    dTdt = -(dqx + dqz);
end

% Function to calculate advection rate
function dTdt = advection2D(f, u0,v0, dx, dz, ind_x, ind_z, ADVN)

% split the velocities into positive and negative
upos_x = 0.5 * (u0 + abs(u0));    % positive velocity
uneg_x = 0.5 * (u0 - abs(u0));    % negative velocity
upos_z = 0.5 * (v0 + abs(v0));    % positive velocity
uneg_z = 0.5 * (v0 - abs(v0));    % negative velocity


% get values on stencil nodes for x-direction
fmmm_x = f(ind_x(1:end-6), :);    % x-stencil
fmm_x  = f(ind_x(2:end-5), :);
fm_x   = f(ind_x(3:end-4), :);
fc_x   = f(ind_x(4:end-3), :);
fp_x   = f(ind_x(5:end-2), :);
fpp_x  = f(ind_x(6:end-1), :);
fppp_x = f(ind_x(7:end), :);

% get values on stencil nodes for z-direction
fmmm_z = f(:, ind_z(1:end-6));    % z-stencil
fmm_z  = f(:, ind_z(2:end-5));
fm_z   = f(:, ind_z(3:end-4));
fc_z   = f(:, ind_z(4:end-3));
fp_z   = f(:, ind_z(5:end-2));
fpp_z  = f(:, ind_z(6:end-1));
fppp_z = f(:, ind_z(7:end));

% calculate heat flux by advection in both directions
switch ADVN
    case 'UPW1'
        fppos_x = fc_x;  fpneg_x = fp_x;
        fmpos_x = fm_x;  fmneg_x = fc_x;

        fppos_z = fc_z;  fpneg_z = fp_z;
        fmpos_z = fm_z;  fmneg_z = fc_z;

    case 'CFD2'
        fppos_x = (fc_x + fp_x) ./ 2;  fpneg_x = fppos_x;
        fmpos_x = (fc_x + fm_x) ./ 2;  fmneg_x = fmpos_x;

        fppos_z = (fc_z + fp_z) ./ 2;  fpneg_z = fppos_z;
        fmpos_z = (fc_z + fm_z) ./ 2;  fmneg_z = fmpos_z;

    case 'UPW3'
        fppos_x = (2*fp_x + 5*fc_x - fm_x) ./ 6;  fpneg_x = (2*fc_x + 5*fp_x - fpp_x) ./ 6;
        fmpos_x = (2*fc_x + 5*fm_x - fmm_x) ./ 6;  fmneg_x = (2*fm_x + 5*fc_x - fp_x) ./ 6;

        fppos_z = (2*fp_z + 5*fc_z - fm_z) ./ 6;  fpneg_z = (2*fc_z + 5*fp_z - fpp_z) ./ 6;
        fmpos_z = (2*fc_z + 5*fm_z - fmm_z) ./ 6;  fmneg_z = (2*fm_z + 5*fc_z - fp_z) ./ 6;

    case 'WENO5'
        fppos_x = weno5poly(fmm_x, fm_x, fc_x, fp_x, fpp_x);  
        fpneg_x = weno5poly(fppp_x, fpp_x, fp_x, fc_x, fm_x); 
        fmpos_x = weno5poly(fmmm_x, fmm_x, fm_x, fc_x, fp_x); 
        fmneg_x = weno5poly(fpp_x, fp_x, fc_x, fm_x, fmm_x);

        fppos_z = weno5poly(fmm_z, fm_z, fc_z, fp_z, fpp_z);  
        fpneg_z = weno5poly(fppp_z, fpp_z, fp_z, fc_z, fm_z); 
        fmpos_z = weno5poly(fmmm_z, fmm_z, fm_z, fc_z, fp_z); 
        fmneg_z = weno5poly(fpp_z, fp_z, fc_z, fm_z, fmm_z);
end

% calculate flux balance for rate of change
div_qpos_x = upos_x .* (fppos_x - fmpos_x) / dx;  % x-direction positive flux
div_qneg_x = uneg_x .* (fpneg_x - fmneg_x) / dx;  % x-direction negative flux
div_q_x    = div_qpos_x + div_qneg_x;           % total x-direction flux

div_qpos_z = upos_z .* (fppos_z - fmpos_z) / dz;  % z-direction positive flux
div_qneg_z = uneg_z .* (fpneg_z - fmneg_z) / dz;  % z-direction negative flux
div_q_z    = div_qpos_z + div_qneg_z;           % total z-direction flux

% align dimensions for combining flux contributions
if size(div_q_x, 1) < size(div_q_z, 1)
    div_q_x = padarray(div_q_x, [size(div_q_z, 1) - size(div_q_x, 1), 0], 'post');
elseif size(div_q_x, 2) < size(div_q_z, 2)
    div_q_z = padarray(div_q_z, [0, size(div_q_x, 2) - size(div_q_z, 2)], 'post');
end

% combine x and z contributions
div_q = div_q_x + div_q_z;

% rate of change
dTdt = -div_q;

end

function [fface] = weno5poly (fmm, fm, fc, fp, fpp)
% 5th order WENO polynomials from Jiang & Shu, 1996, J Comp Physics

% 5th order polynomials
p1 = (2*fmm - 7*fm + 11*fc )/6;
p2 = ( -fm  + 5*fc +  2*fp )/6;
p3 = (2*fc  + 5*fp -    fpp)/6;

% smoothness measure
b1 = 13/12*(fmm - 2*fm + fc ).^2 + 1/4*(  fmm - 4*fm + 3*fc ).^2;
b2 = 13/12*(fm  - 2*fc + fp ).^2 + 1/4*(  fm  -          fp ).^2;
b3 = 13/12*(fc  - 2*fp + fpp).^2 + 1/4*(3*fc  - 4*fp +   fpp).^2;

% weights
g   = [1/10, 6/10, 3/10];
wp1 = g(1)./(b1.^2 + eps);
wp2 = g(2)./(b2.^2 + eps);
wp3 = g(3)./(b3.^2 + eps);

% get flux (normalize weights at the same time)
fface = (wp1.*p1 + wp2.*p2 + wp3.*p3) ./ (wp1 + wp2 + wp3) ;

end