% Description:
%   Compute OAD/OED design with fixed budget C using rounding & extension

%% Initialization
criterion = 'A';  % {'D', 'A', 'c', 'Ds'}
C = 100;            % Total cost budget
M = 150;
q_cost = 0.2;
addpath('../util');

cVec_c  = [0; 1; 1];
cVec_Ds = [1; 0; 0];
cVec_struct.cVec_c  = cVec_c;
cVec_struct.cVec_Ds = cVec_Ds;
S = [1, M]';      % Design space
theta = [0.07; 0.93; 0.96];
tol = 1E-5;
q = length(theta);
u = S(1):S(2);
N = length(u);
zero_vec = zeros(N,1);
one_vec  = ones(N,1);

%% 2. CVX: Approximate Optimal Design
cvx_begin quiet
cvx_precision best
variables w(N,1) del(1)
minimize del
subject to
    M = compute_FIM_GT_cost(u, w, theta, q_cost);
    switch criterion
        case 'D'
            -log_det(M) <= del;
        case 'A'
            trace_inv(M) <= del;
        case 'c'
            matrix_frac(cVec_c, M) <= del;
        case 'Ds'
            matrix_frac(cVec_Ds, M) <= del;
        otherwise
            error("Invalid criterion");
    end
    w >= 0;
    one_vec' * w == 1;
cvx_end

kk = find(w > tol);
design_app = [u(kk); w(kk)'];
xx = design_app(1,:);
ww = design_app(2,:);
cxx = 1 - q_cost + q_cost * xx;

%% 3. Floor allocation under budget C
nn_i = C * ww ./ cxx;
n_prime = floor(nn_i);
design_round1 = [xx; n_prime];
Cr = round(C - sum(n_prime .* cxx), 4);  % remaining cost

%% 4. Extended Support Points
n_index = 2;
support_ext = [];

for i = 1:length(xx)
    support_ext = [support_ext, (xx(i) - n_index):(xx(i) + n_index)];
end
xx_temp = unique(support_ext);
xx_temp = xx_temp(ismember(xx_temp, u));
cxx_all = 1 - q_cost + q_cost * xx_temp;

%% 5. All Valid Rounding Extensions
des_exact_rounding2 = calc_budget_round_combinations(cxx_all, Cr);

%% 6. Combine Rounded + Extended Designs
x = design_round1;
n_cases = size(des_exact_rounding2, 1);
result_cell = cell(n_cases, 1);

for k = 1:n_cases
    x_ext = xx_temp;
    w_ext = des_exact_rounding2(k, 1:end-2);
    all_keys = union(x(1,:), x_ext);
    values = zeros(1, length(all_keys));

    [~, loc_x] = ismember(x(1,:), all_keys);
    values(loc_x) = values(loc_x) + x(2,:);

    [~, loc_y] = ismember(x_ext, all_keys);
    values(loc_y) = values(loc_y) + w_ext;

    % Filter out 0 weights
    nz_idx = values > 0;
    result_cell{k} = [all_keys(nz_idx); values(nz_idx)];
end

%% 7. Evaluate All Extended Designs
loss_round = zeros(n_cases, 1);

for i = 1:n_cases
    des = result_cell{i};
    cweights = 1 - q_cost + q_cost * des(1,:);
    w_norm = des(2,:) .* cweights / C;
    M = compute_FIM_GT_cost(des(1,:), w_norm, theta, q_cost);

    switch criterion
        case 'D'
            loss_round(i) = calc_loss_D(M);
        case 'A'
            loss_round(i) = calc_loss_A(M);
        case 'c'
            loss_round(i) = calc_loss_c(M, cVec_c);
        case 'Ds'
            loss_round(i) = calc_loss_c(M, cVec_Ds);
    end
end

%% 8. Rank and Select Best Design
des_exact_rounding2(:, end+1) = loss_round;
sorted_des = sortrows(des_exact_rounding2, size(des_exact_rounding2, 2));

x_names = arrayfun(@(v) sprintf('%.0f', v), xx_temp, 'UniformOutput', false);
column_names = [x_names, {'Used Cost', 'Remaining', 'Loss'}];
T = array2table(sorted_des, 'VariableNames', column_names);

disp('Rounded Floor Allocation:');
disp(design_round1);

disp('Top 5 Rounded Extended Designs:');
disp(T(1:5, :));

%% 9. Final Design: Combine Floor + Best Extended Design
best_ext = table2array(T(1, 1:end-3));
y = [xx_temp; best_ext];
all_keys = union(x(1,:), y(1,:));
combined_values = zeros(size(all_keys));

[~, loc_x] = ismember(x(1,:), all_keys);
[~, loc_y] = ismember(y(1,:), all_keys);

combined_values(loc_x) = combined_values(loc_x) + x(2,:);
combined_values(loc_y) = combined_values(loc_y) + y(2,:);
nz_idx = combined_values > 0;

design_filter = [all_keys(nz_idx); combined_values(nz_idx)];
cxx_filter = 1 - q_cost + q_cost * design_filter(1,:);
w_final = design_filter(2,:) .* cxx_filter / C;

%% 10. Efficiency Evaluation
M_ex = compute_FIM_GT_cost(design_filter(1,:), w_final, theta, q_cost);
M_app = compute_FIM_GT_cost(design_app(1,:), design_app(2,:), theta, q_cost);

switch criterion
    case 'D'
        loss_ex = calc_loss_D(M_ex);
        loss_app = calc_loss_D(M_app);
        eff = (loss_app/loss_ex)^(1/q);
    case 'A'
        loss_ex = calc_loss_A(M_ex);
        loss_app = calc_loss_A(M_app);
        eff =  loss_app/ loss_ex;
    case 'c'
        loss_ex = calc_loss_c(M_ex, cVec_c);
        loss_app = calc_loss_c(M_app, cVec_c);
        eff = loss_app/loss_ex;
    case 'Ds'
        loss_ex = calc_loss_c(M_ex, cVec_Ds);
        loss_app = calc_loss_c(M_app, cVec_Ds);
        eff = loss_app/loss_ex;
end

%% 11. Final Output
total_cost_used = sum(design_filter(2,:) .* (1 - q_cost + q_cost * design_filter(1,:)));

fprintf('\n--- Final Budgeted Design Summary ---\n');
disp('Final Rounded Design (Exact Runs):');
disp(design_filter);

fprintf('Efficiency (%s-optimality): %.4f\n', criterion, eff);
fprintf('Total Cost Used: %.4f / %.2f | Remaining: %.4f\n', ...
    total_cost_used, C, C - total_cost_used);