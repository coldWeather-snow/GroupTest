%Compute the maximin design with dual objective for group testing
%experiment

%% addpath
% clear;
start_time = tic;  % Start timer
addpath('util');

%% Set up the parameters
q_cost = 0.2;
tol = 1E-5;
M = 61;
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
% result_Ds= compute_design_SO(u, theta, q_cost, q, 'Ds', cVec_struct, tol);
result_c = compute_design_SO(u, theta, q_cost,  q, 'c', cVec_struct, tol);
result_E = compute_design_SO(u, theta, q_cost,  q, 'E', cVec_struct, tol);
% loss_single =  struct('D', result_D.loss, 'A', result_A.loss, 'Ds', result_Ds.loss, 'c', result_c.loss,...
  % 'E', result_E.loss);
loss_single =  struct('D', result_D.loss, 'A', result_A.loss, 'c', result_c.loss,...
  'E', result_E.loss);

%%% Four criterion
% q_cost = 0;
result_DEDsc = compute_maximin_design2(u, theta, q_cost,  q, loss_single, {'D', 'E', 'A', 'c'}, ...
  cVec_struct, tol);
round(result_DEDsc.design, 3)
result_DEDsc.efficiency
round(1/result_DEDsc.tstar,3)

elapsed_time = toc(start_time);
fprintf('Elapsed time: %.3f seconds\n', elapsed_time);


result_EA = compute_maximin_design2(u, theta, q_cost,  q, loss_single, {'E', 'A'}, cVec_struct, tol);

result_ED = compute_maximin_design2(u, theta, q_cost,  q, loss_single, {'E', 'D'}, cVec_struct, tol);
result_Ec = compute_maximin_design2(u, theta, q_cost,  q, loss_single, {'E', 'c'}, cVec_struct, tol);
result_EDs = compute_maximin_design2(u, theta, q_cost,  q, loss_single, {'E', 'Ds'}, cVec_struct, tol);
1/result_EA.tstar
1/result_ED.tstar
1/result_Ec.tstar
1/result_EDs.tstar

%% calculate derivative 

dd_ED = calc_directional_derivatives(u, result_ED.M, theta, q_cost, {'E', 'D'});
eta_ED = calc_eta_weights2(result_ED.tstar, loss_single, result_ED.loss, dd_ED, {'E', 'D'}, tol);

%%% Display the result
fprintf("The value of eta's are: ");
fprintf('%.4f ', eta_ED);  % adjust number format as needed
fprintf('\n');

fprintf("The value of tstar is: ");
fprintf('%.4f ', result_ED.tstar);  % adjust number format as needed
fprintf('\n');

%%% compute the directional derivative for multi-objectivity.
d_multi = eta_ED(1) * dd_ED.dE + eta_ED(2) * dd_ED.dD;
figure;
plot(d_multi)
hold on;

% Add vertical lines at the design points
x_vals = result_ED.design(1, :);
for i = 1:length(x_vals)
    xline(x_vals(i), '--r', 'LineWidth', 1.5);
end

% Add a vertical line at y = 0
yline(0, '--k', 'LineWidth', 1.5); % Dashed black line
hold off;


dd_ED = calc_directional_derivatives(u, result_ED.M, theta, q_cost, {'E', 'D'});
eta_ED = calc_eta_weights2(result_ED.tstar, loss_single, result_ED.loss, dd_ED, {'E', 'D'}, tol);

%%% Display the result
fprintf("The value of eta's are: ");
fprintf('%.4f ', eta_ED);  % adjust number format as needed
fprintf('\n');

fprintf("The value of tstar is: ");
fprintf('%.4f ', result_ED.tstar);  % adjust number format as needed
fprintf('\n');

%%% compute the directional derivative for multi-objectivity.
d_multi = eta_ED(1) * dd_ED.dE + eta_ED(2) * dd_ED.dD;
figure;
plot(d_multi)
hold on;

% Add vertical lines at the design points
x_vals = result_ED.design(1, :);
for i = 1:length(x_vals)
    xline(x_vals(i), '--r', 'LineWidth', 1.5);
end

% Add a vertical line at y = 0
yline(0, '--k', 'LineWidth', 1.5); % Dashed black line
hold off;


%% EA
dd_EA = calc_directional_derivatives(u, result_EA.M, theta, q_cost, {'E', 'A'});
eta_EA = calc_eta_weights2(result_EA.tstar, loss_single, result_EA.loss, dd_EA, {'E', 'A'}, tol);

%%% Display the result
fprintf("The value of eta's are: ");
fprintf('%.4f ', eta_EA);  % adjust number format as needed
fprintf('\n');

fprintf("The value of tstar is: ");
fprintf('%.4f ', result_EA.tstar);  % adjust number format as needed
fprintf('\n');

%%% compute the directional derivative for multi-objectivity.
d_multi_EA = eta_EA(1) * dd_EA.dE + eta_EA(2) * dd_EA.dA;
figure;
plot(d_multi_EA)
hold on;

% Add vertical lines at the design points
x_vals = result_EA.design(1, :);
for i = 1:length(x_vals)
    xline(x_vals(i), '--r', 'LineWidth', 1.5);
end

% Add a vertical line at y = 0
yline(0, '--k', 'LineWidth', 1.5); % Dashed black line
hold off;

