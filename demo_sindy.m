% DEMO_SINDY Example calls for the function version of SINDy/SSR.
%
% Run this file from the root folder of the repository.

clear; clc; close all;

file = 'trajectories/limit_cycle_lambda_-15.000.mat';

% Learn the gradient component: F = -p/2.
results_gradient = sindy_fct(file, 'gradient', true);


% Learn the transverse component: F = dotX - p/2.
results_transverse = sindy_fct(file, 'transverse', true);
