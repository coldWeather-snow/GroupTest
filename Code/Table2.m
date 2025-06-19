%Compute the maximin design with dual objective for group testing
%experiment

%% addpath
% clear;
start_time = tic;  % Start timer
addpath('util');

%% Set up the parameters
q_cost = 0;
tol = 1E-5;
M = 150;
S = [1, M]'; % design space
cVec_c = [0, 1, 1]';
cVec_Ds = [1, 0, 0]';
cVec_struct = struct('cVec_c', cVec_c, 'cVec_Ds', cVec_Ds);
theta = [0.07, 0.93, 0.96]';
p0 = theta(1); p1 = theta(2); p2 = theta(3);
q = length(theta);
u =  S(1) : S(2);
N = length(u);

%% Compute the single objectives for GT (group testing) 
result_D = compute_design_SO(u, theta, q_cost, q, 'D', cVec_struct, tol);
result_A = compute_design_SO(u, theta, q_cost, q, 'A', cVec_struct, tol);
result_Ds= compute_design_SO(u, theta, q_cost, q, 'Ds', cVec_struct, tol);
result_c = compute_design_SO(u, theta, q_cost,  q, 'c', cVec_struct, tol);
result_E = compute_design_SO(u, theta, q_cost,  q, 'E', cVec_struct, tol);
loss_single =  struct('D', result_D.loss, 'A', result_A.loss, 'Ds', result_Ds.loss, 'c', result_c.loss,...
  'E', result_E.loss);

result_DA = compute_maximin_design(u, theta, q_cost,  q, loss_single, {'D', 'A'}, cVec_struct, tol);
result_DDs = compute_maximin_design(u, theta, q_cost, q, loss_single, {'D', 'Ds'}, cVec_struct, tol);
result_DADs = compute_maximin_design(u, theta, q_cost, q, loss_single, {'D', 'A', 'Ds'}, cVec_struct, tol);
result_DDsc = compute_maximin_design(u, theta, q_cost, q, loss_single, {'D', 'Ds', 'c'}, cVec_struct, tol);

%% Output the results

%%% Single-criterion
result_D.design
result_A.design
result_Ds.design
result_c.design

%%% dual-criterion
result_DA.design
result_DDs.design

%%% Triple-criterion
result_DADs.design
result_DDsc.design

%%% Display the 1/t_star, note this is only available in the
%%% multi-objective design
fprintf('q_cost = %.1f with M = %d\n', q_cost, M);
fprintf('1/tstar for DA   = %.4f\n', 1 / result_DA.tstar);
fprintf('1/tstar for DDs  = %.4f\n', 1 / result_DDs.tstar);
fprintf('1/tstar for DADs = %.4f\n', 1 / result_DADs.tstar);
fprintf('1/tstar for DDsc = %.4f\n', 1 / result_DDsc.tstar);
