%Session 4
%One-dimensional (1D) model of heat transport by advection and diffusion.
%Implicit Time Integration,
%***** RUN 1D ADVECTION DIFFUSION MODEL ***********************************

% clear workspace
clear all; close all; %clc;

%NN = [100,200,400];

%for nn = 1:3

% set model parameters
W  = 500;          % domain width [m]
w  = 50;           % magma width [m]
N  = 200           % grid size
dx = W/N;          % grid spacing
%N  = NN(nn);      % grid size (for loop through)

T0   = 100;        % initial background temperature [C]
dT   = 1000;       % initial temperature peak amplitude [C]
wT   = 10;         % initial temperature peak width [m]
k0   = 1e-6;       % heat diffusivity [m2/s]
u0   = 1e-6;       % advection speed [m/s]
BC   = 'periodic'; % boundary condition option flag ('insulating', 'periodic')
ADVN = 'WENO5';    % advection scheme ('UPW1', 'CFD2', 'UPW3', 'WENO5')
TINT = 'BE1';      % time integration scheme ('FE1','RK2','HE2','RK4','BE1','CN2')

yr    = 3600*24*365;  % seconds per year [s]
tend  = W/u0;         % stopping time [s]
CFL   = 1/100;        % Time step limiter
nop   = 1000;         % output figure produced every 'nop' steps

%*****  RUN MODEL
run('./advdiff1D_implicit.m');

%end