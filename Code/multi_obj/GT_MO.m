%Compute the maximin design with dual objective for group testing
%experiment
clear;
addpath('util');
% addpath

q_cost = 0.8;
tol = 1E-5;
S = [1, 150]'; % design space
cVec_c = [0, 1, 1]';
cVec_Ds = [1, 0, 0]';
cVec_struct = struct('cVec_c', cVec_c, 'cVec_Ds', cVec_Ds);
theta = [0.07, 0.93, 0.96]';
p0 = theta(1); p1 = theta(2); p2 = theta(3);
q = length(theta);
u =  S(1) : S(2);
N = length(u);

%% Compute the single objectives for GT 
result_D = compute_design_SO(u, theta, q_cost, q, 'D', cVec_struct, tol);
result_A = compute_design_SO(u, theta, q_cost, q, 'A', cVec_struct, tol);
result_Ds= compute_design_SO(u, theta, q_cost, q, 'Ds', cVec_struct, tol);
result_c = compute_design_SO(u, theta, q_cost,  q, 'c', cVec_struct, tol);
loss_single =  struct('D', result_D.loss, 'A', result_A.loss, 'Ds', result_Ds.loss, 'c', result_c.loss);

result_DA = compute_maximin_design(u, theta, q_cost,  q, loss_single, {'D', 'A'}, cVec_struct, tol);
result_DDs = compute_maximin_design(u, theta, q_cost, q, loss_single, {'D', 'Ds'}, cVec_struct, tol);
result_DADs = compute_maximin_design(u, theta, q_cost, q, loss_single, {'D', 'A', 'Ds'}, cVec_struct, tol);
result_DDsc = compute_maximin_design(u, theta, q_cost, q, loss_single, {'D', 'Ds', 'c'}, cVec_struct, tol);

%%% Single
result_D.design
result_A.design
result_Ds.design
result_c.design

%%% daul 
result_DA.design
result_DDs.design

%%% Triple
result_DADs.design
result_DDsc.design

%% Equivalence theorem 
% D- & A-
dd_DA = calc_directional_derivatives(u, result_DA.M, theta, q_cost, {'D', 'A'});
eta_DA = calc_eta_weights(result_DA.tstar, loss_single, result_DA.loss, dd_DA, {'D', 'A'}, tol);

%%% compute the directional derivative for multi-objectivity.
d_multi = eta_DA(1) * dd_DA.dD + eta_DA(2) * dd_DA.dA;
figure;

% --- First plot: D-optimality ---
subplot(1, 3, 1);
plot(dd_DA.dD, 'b', 'LineWidth', 2)
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
title('(a)', 'FontSize', 14);
xlabel('$$u_i$$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$$d_D(u_i, w^{**})$$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;

% --- Second plot: A-optimality ---
subplot(1, 3, 2);
plot(dd_DA.dA, 'b', 'LineWidth', 2)
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
title('(b)', 'FontSize', 14);
xlabel('$$u_i$$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$$d_A(u_i, w^{**})$$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;

% --- Third plot: Maximin (D+A) ---
subplot(1, 3, 3);
plot(d_multi, 'b', 'LineWidth', 2)
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
title('(c)', 'FontSize', 14);
xlabel('$$u_i$$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$$\sum_{j=1}^2 \eta_j d_j(u_i, \mathbf{w}^{**})$$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;


%%%  D-, A- and Ds
dd_DADs = calc_directional_derivatives(u, result_DADs.M, theta, q_cost, {'D', 'A', 'Ds'});
eta_DADs = calc_eta_weights(result_DADs.tstar, loss_single, result_DADs.loss, dd_DADs, {'D', 'A', 'Ds'}, 1E-5);

%%% compute the directional derivative for multi-objectivity.
d_multi_DADs = eta_DADs(1) * dd_DADs.dD + eta_DADs(2) * dd_DADs.dA + eta_DADs(3) *dd_DADs.dDs;
figure;

% --- First plot: D-optimality (DADs) ---
subplot(2, 2, 1);
plot(dd_DADs.dD, 'b', 'LineWidth', 2)
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
title('(a)', 'FontSize', 14);
xlabel('$$u_i$$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$$d_D(u_i, \mathbf{w}^{**})$$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;

% --- Second plot: A-optimality (DADs) ---
subplot(2, 2, 2);
plot(dd_DADs.dA, 'b', 'LineWidth', 2)
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
title('(b)', 'FontSize', 14);
xlabel('$$u_i$$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$$d_A(u_i, \mathbf{w}^{**})$$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;

% --- Third plot: Ds-optimality (DADs) ---
subplot(2, 2, 3);
plot(dd_DADs.dDs, 'b', 'LineWidth', 2)
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
title('(c)', 'FontSize', 14);
xlabel('$$u_i$$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$$d_{D_s}(u_i, \mathbf{w}^{**})$$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;

% --- Fourth plot: Maximin (D+A) ---
subplot(2, 2, 4);
plot(d_multi_DADs, 'b', 'LineWidth', 2)
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
title('(d)', 'FontSize', 14);
xlabel('$$u_i$$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$$\sum_{j=1}^3 \eta_j d_j(u_i, \mathbf{w}^{**})$$', 'Interpreter', 'latex', 'FontSize', 14);
grid on;

eta_DADs