%Session 1
%Making a very simple thermal evolution model of the Earths Interior 

%The mantle of the young Earth was initaly much hotter than it is today. It
%has been loosing heat to space throguh the atmosphere, some of which is
%being repenished by heat drawn from the even hotter core as well as from
%radioactive decay of some unstable elements in its composition. Lets build
%a very simple model to simulate the thermal evolution of the Earths mantle
%through time. We start by casting the conceptual model decribed above inta

% differential equation

%dT_m / dt = q_ma + qcm + Q_r

%Where
%T_m :  mantle temperatuer
%d/dt :  time derivative,
%r_ma : flux of heat from mantle to the atmosphere
%q_cm : flux of heat from the core to the mantle
%Q_r : radioactive heating rate.

%model the heat flow rates from the core to the mantle and the mantle to
%the atmosphere as proportional to their temperature difference and modulated
%by an effective heat transport coefficient . 

% Qma = K(Ta - Tm)
% Qmc = K(Tc - Tm)


%To keep things simple we can say that heat transport within the liquid outer
%core and gaseous atmosphere are much more rapid than in the solid state mantle
%and hence the rate of heat transfer will be limited by the overturn rate of
%mantle convection. We therefore relate the rate coefficent the characteristic
%time of mantle overturn.

%%
clear; close all; clc; % clear workspace, close variales, clean command window

%Define Variables 
Tm = 2500;          % the inital temperature of the mantle (K)
Tc = 3000;          % the constant temperature of the core (K)
Ta = 300;           % the constant atmosphere temperature (K)
Qr = 6e-14;         % the inital radioactive rate (K/s)

yr = 365*24*3600;   %seconds per year [s]
tend = 4.5e9 *yr;   %sets a stopping time
t = 0;              %sets inital time 
tau = 5e8 * yr;     %the mantle overturn time [s]
K = 1/tau;          %heat transfer coefficent [1/s]
dt = 5e7 * yr;      %defines the time step [s]


%loop through time untill the stopping time is reached
while t <= tend

    %calculate the heat transfers

    qMa = K*(Ta-Tm);
    qMc = K*(Tc-Tm);
    %update the mantle temperature
    Tm = Tm + (qMa + qMc + Qr) * dt; %This calculates the Tm after a time step
    
    %update the core temperature
    Tc = Tc -qMc *dt;

    %Displays the time, the temperatuer of the mantle and the temperature
    %of the core.
    disp(['time = ' num2str(t/1e6/yr) ' Myr   Tm = ' num2str(Tm) 'Tc = ' num2str(Tc) ' K'])


    %plot model progress
    plot(t/yr/1e6,Tm,'ro',t/yr/1e6,Tc,'bo'); axis tight; hold on;drawnow;
    legend('mantle T', 'core T')
    xlabel('Time [Myr]');
    ylabel('Temperature [K]');

     %increment time 
    t = t + dt;

end 
