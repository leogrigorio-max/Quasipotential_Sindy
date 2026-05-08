% ==================================================
% Initializes variables and calls functions
% ==================================================

clear
clc
addpath NR_functions
N = 4e2;

s = linspace(0,1,N); ds = s(2)-s(1);

epsilon = 0.25;

delta = 1.0; %initial relative error

x = zeros(2,N);
p = zeros(2,N);

x_ = zeros(2,N);
p_ = zeros(2,N);

x_action = [];


C = 1.0; %initial arclength guess
