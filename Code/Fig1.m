%Compute the maximin design with dual objective for group testing
%experiment

% addpath
clear;
addpath('util');
run('Table2.m') 
% we obtain the results of the multi-objective group testing designs from
% this script, make sure to change the value of the q_cost variable and M

%% Equivalence theorem 
% D- & A-
dd_DA = calc_directional_derivatives(u, result_DA.M, theta, q_cost, {'D', 'A'});
eta_DA = calc_eta_weights(result_DA.tstar, loss_single, result_DA.loss, dd_DA, {'D', 'A'}, tol);

%%% Display the result
fprintf("The value of eta's are: ");
fprintf('%.4f ', eta_DA);  % adjust number format as needed
fprintf('\n');

fprintf("The value of tstar is: ");
fprintf('%.4f ', result_DA.tstar);  % adjust number format as needed
fprintf('\n');

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
