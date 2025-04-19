%Compute the maximin design with dual objective for group testing
%experiment

% addpath
clear;
addpath('util');
run('Table2.m') 
% we obtain the results of the multi-objective group testing designs from
% this script, make sure to change the value of the q_cost variable and M

%% Equivalence theorem 
%%%  D-, A- and Ds
dd_DADs = calc_directional_derivatives(u, result_DADs.M, theta, q_cost, {'D', 'A', 'Ds'});
eta_DADs = calc_eta_weights(result_DADs.tstar, loss_single, result_DADs.loss, dd_DADs, {'D', 'A', 'Ds'}, 1E-5);

%%% Display the result
fprintf("The value of q_cost is: ");
fprintf('%.1f ', q_cost);  
fprintf('\n');
fprintf("The value of M is: ");
fprintf('%.0f ', M); 
fprintf('\n');

fprintf("The value of eta's are: ");
fprintf('%.4f ', eta_DADs); 
fprintf('\n');

fprintf("The value of tstar is: ");
fprintf('%.4f ', result_DADs.tstar);  
fprintf('\n');



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