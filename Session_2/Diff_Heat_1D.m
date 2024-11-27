%Session 2
%One-dimensional (1D) model of heat transport by diffusion.
%Representing a cooling magma filled fracture within a cooler wall rock.


%%
clear; close all; clc; % clear workspace, close variales, clean command window

%set model parameters 

W  = 500;     % domain width [m]
w  = 50;      % magma width [m]
N  = 200;     % grid size
dx = W/N;     % grid spacing

Tr = 100;     % initial rock temperature [C]
Tm = 1100;    % initial magma temperature [C]
kr = 1e-6;    % rock diffusivity [m2/s]
km = 1e-5;    % magma diffusivity [m2/s]
BC = 'insulating'; % boundary condition option flag ('isothermal' or 'insulating')

yr    = 3600*24*365;      % seconds per year [s]
tend  = 100*yr;           % stopping time [s]
dt    = (dx/2)^2/max([kr,km]);  % time step [s]
t     = 0;                % initial time [s]
k     = 0;                % initial time step count

% create x-coordinate vectors
xc = dx/2:dx:W-dx/2;    % coordinate vector for cell centre positions [m]
xf = 0:dx:W;            % coordinate vectore for cell face positions [m]

% set variable diffusion coefficients at cell centres
kv     = zeros(size(xc)) + kr;                 % initialise variable k array
ind    = find(xc >= W/2-w/2 & xc <= W/2+w/2);  % find indeces for magma domain
kv(ind)= km;                                   % set magma domain
kv     = [kr, kv, kr];                         % add ghost nodes
kf     = (kv(1:end-1) + kv(2:end))/2;          % get mid-point values at cell faces

% set initial condition for temperature at cell centres
T      = zeros(size(xc)) + Tr;                 % initialise T array at Tr
T(ind) = Tm;                                   % set magma domain to Tm

% add ghost nodes for boundary conditions
T  = [Tr, T, Tr];    % add ghost nodes at both boundaries
ic = 2:N+1;          % indeces of T points interior to domain
T0 = T;              % store initial condition

% initialise figure for results
figure(); clf

while t <= tend
    % calculate heat flux
    q = - kf .* diff(T)/dx;

    % calculate rate of change
    dTdt = - diff(q)/dx;

    % update temperature (interior points only)
    T(ic) = T(ic) + dTdt * dt;

    % enforce boundary conditions
    switch BC
        case 'isothermal'
            T([1 end]) = Tr;
        case 'insulating'
            T([1 end]) = T([2 end-1]);
    end

    % plot model progress
    if ~mod(k,10)
    plot(xc,T(ic),'r-',xc,T0(ic),'k-','LineWidth',1.5); axis tight; box on;
    xlabel('x [m]','FontSize',15)
    ylabel('T [C]','FontSize',15)
    title(['Evolving Temperature; time = ',num2str(t/yr),' yr'],'FontSize',18)
    drawnow;
    end

    % increment time and step count
    t = t+dt;
    k = k+1;
end


%%
%Heat Flux Calculation Versions:

%%% Version 1: For loop %%%

%A for loop can be used to calcuate the fluxes at the cell faces and store
%the in a vector q of length N+1.

% q = zeros(1,N+1);
% for i=1:N+1
%     q(i) = - k0 * (T(i+1) - T(i))/dx;
% end

%%
%%% Version 2: Vectorisation %%%


%The loops in matlab can be inefficent and therfore its more efficient to
%vectorise a calculation. This involves commands which perform operations on 
% whole arrays at once. THis can be achieved using element-wise (eg .*) or
% matrix (*) operations. 

%q = - k0 * (T(2:end) - T(1:end-1))/dx;

%%
%%% Version 3: diff() %%%


%The final version is to use the diff()function. Also acounting for the
%time dependent model.

% calculate heat flux
%q = - k0 * diff(T)/dx;
% calculate temperature rate of change
%dTdt = - diff(q)/dx;

%OR 

%This can be written as one line like so

% dTdt = k0 * diff(T,2)/dx^2;
