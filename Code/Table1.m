% This script is used to reproduce the results given in Table 1
% it is about single objective design 
% with M = 61, with various q_cost \in [0,1]

%% Add utility functions from util directory
addpath('../util');

%% Initialzie the paramters 
M = 61;
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

% respective loss functions
disp("Loss function for each criterion at the optimal designs")
loss_Vec = [res_D.loss, res_A.loss, res_Ds.loss, res_c.loss]';
T_loss = array2table(loss_Vec', 'VariableNames', {'D-', 'A-', 'Ds-', 'c-'});
T_loss.Properties.RowNames = {'loss_opt'};
disp(T_loss);

% calculate the relative efficiency
disp("Efficiency of D-criterion");
loss_DonD = calc_loss_D(res_D.M);
loss_DonA = calc_loss_A(res_D.M);
loss_DonDs = calc_loss_c(res_D.M, cVec_Ds);
loss_Donc = calc_loss_c(res_D.M, cVec_c);
eff_D_criterion = loss_Vec' ./ [loss_DonD, loss_DonA, loss_DonDs, loss_Donc];
T_eff = array2table(eff_D_criterion, 'VariableNames', {'Eff_D', 'Eff_A', 'Eff_Ds', 'Eff_c'});
T_eff.Properties.RowNames = {'D-criterion'};
disp(T_eff);