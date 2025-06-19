% This script is used to reproduce the results given in Table 1
% it is about single objective design 
% with M = 61, with various q_cost \in [0,1]

%% Add utility functions from util directory
addpath('./util');
%% Initialzie the paramters 
M = 150;
S = [1, M]'; % design space
 theta = [0.07, 0.93, 0.96]';
q_cost = 0.8;
cVec_c = [0, 1, 1]';
cVec_Ds = [1, 0, 0]';
p0 = theta(1); p1 = theta(2); p2 = theta(3);
q = length(theta);
u =  S(1) : S(2);
N = length(u);
tol = 1E-5;
cVec_struct = struct('cVec_c', cVec_c, 'cVec_Ds', cVec_Ds);

%% Compute single objective optimal approximate designs
disp("D-optimality")
res_D = compute_design_SO(u, theta, q_cost, q, 'D', cVec_struct, tol);
disp("A-optimality")
res_A = compute_design_SO(u, theta, q_cost, q, 'A', cVec_struct, tol);
disp("Ds-optimality")
res_Ds = compute_design_SO(u, theta, q_cost, q, 'Ds', cVec_struct, tol);
disp("c-optimality")
res_c = compute_design_SO(u, theta, q_cost, q, 'c', cVec_struct, tol);
disp("E-optimality")
res_E = compute_design_SO(u, theta, q_cost, q, 'E', cVec_struct, tol);

%% Analying the results
% resulting optimal approximate designs
disp("D-optimality")
res_D.design
disp("A-optimality")
res_A.design
disp("Ds-optimality")
res_Ds.design
disp("c-optimality")
res_c.design
disp("E-optimality")
res_E.design

% respective loss functions
disp("Loss function for each criterion at the optimal designs")
loss_Vec = round([res_D.loss, res_A.loss, res_Ds.loss, res_c.loss, res_E.loss], 3)';
T_loss = array2table(loss_Vec', 'VariableNames', {'D-', 'A-', 'Ds-', 'c-', 'E'});
T_loss.Properties.RowNames = {'loss_opt'};
disp(T_loss);

% calculate the relative efficiency
calc_eff_D = @(res_opt, res_other) (calc_loss_D(res_opt.M) / calc_loss_D(res_other.M))^(1/3);
calc_eff_E = @(res_opt, res_other) calc_loss_E(res_other.M)/ calc_loss_E(res_opt.M);
calc_eff_Ds = @(res_opt, res_other) calc_loss_c(res_opt.M, cVec_Ds)/calc_loss_c(res_other.M, cVec_Ds);
calc_eff_c = @(res_opt, res_other) calc_loss_c(res_opt.M, cVec_c)/calc_loss_c(res_other.M, cVec_c);
calc_eff_A = @(res_opt, res_other) calc_loss_A(res_opt.M)/calc_loss_A(res_other.M);

res_target = res_E;
eff_onD = calc_eff_D(res_D, res_target);
eff_onE = calc_eff_E(res_E, res_target);
eff_onDs = calc_eff_Ds(res_Ds, res_target);
eff_onA = calc_eff_A(res_A, res_target);
eff_onc = calc_eff_c(res_c, res_target);

eff_all = round([eff_onD, eff_onA, eff_onDs, eff_onc, eff_onE], 3);
T_eff = array2table(eff_all, 'VariableNames', {'Eff_onD', 'Eff_onA', 'Eff_onDs', 'Eff_onc', 'Eff_onE'});
disp(T_eff);