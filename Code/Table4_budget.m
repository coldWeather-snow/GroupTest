% Description:
%   Multi-objective rounding with a fixed total budget (cost C)


%% get result from Multi-objective design
clear;
addpath('util');
run('Table3.m')  % remember to change the M and q_cost in Table3.m

%% 1. Setup
C = 1000;                        % Total cost budget
% q_cost = 0.2;
% theta = [0.07; 0.93; 0.96];

% From prior maximin design
% criteria = {'D', 'A'};
% design_app = result_DA;


% criteria = {'D', 'Ds'};
% design_app = result_DDs;

% criteria = {'D', 'Ds', 'c'};
% design_app = result_DDsc;


criteria = {'D', 'A', 'Ds'};
design_app = result_DADs;

% % Criteria-specific vectors
% cVec_struct.cVec_Ds = [1; 0; 0];
% cVec_struct.cVec_c  = [0; 1; 1];

% Base design
xx = design_app.design(1,:);
ww = design_app.design(2,:);
cxx = 1 - q_cost + q_cost * xx;

%% 2. First-Round (Cost-aware Floor Allocation)
nn_i = C * ww ./ cxx;
n_prime = floor(nn_i);

% Fix zero allocations using redistribution
deficit_idx = find(n_prime == 0);
excess_idx = find(n_prime >= 2);
for i = 1:length(deficit_idx)
    if isempty(excess_idx)
        error('Not enough excess runs to fix zero entries.');
    end
    from = excess_idx(1); to = deficit_idx(i);
    n_prime(from) = n_prime(from) - 1;
    n_prime(to)   = n_prime(to) + 1;
    if n_prime(from) < 2
        excess_idx(1) = [];
    end
end

design_round1 = [xx; n_prime];
Cr = round(C - sum(n_prime .* cxx), 4);  % Remaining budget

%% 3. Extended Support Points
n_index = 2;
extended_pts = [];
for i = 1:length(xx)
    extended_pts = [extended_pts, (xx(i) - n_index):(xx(i) + n_index)];
end
% u = S(1):S(2);  % Assuming S is still in workspace
xx_temp = unique(extended_pts);
xx_temp = xx_temp(ismember(xx_temp, u));
cxx_all = 1 - q_cost + q_cost * xx_temp;

%% 4. Generate Valid Cost-Constrained Extensions
des_exact_rounding2 = calc_budget_round_combinations(cxx_all, Cr);

%% 5. Combine Base + Extensions into Designs
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

    nz_idx = values > 0;
    result_cell{k} = [all_keys(nz_idx); values(nz_idx)];
end

%% 6. Evaluate Designs (Loss + Efficiency)
loss_struct_all = cell(n_cases, 1);
eff_struct_all = cell(n_cases, 1);
eff_scores = zeros(n_cases, 1);

for i = 1:n_cases
    des = result_cell{i};
    des = des(:, des(2,:) > 0);
    cweights = 1 - q_cost + q_cost * des(1,:);
    des(2,:) = des(2,:) .* cweights / C;

    M = compute_FIM_GT_cost(des(1,:), des(2,:), theta, q_cost);
    loss_i = compute_losses(M, criteria, cVec_struct);
    eff_i = compute_efficiencies(loss_single, loss_i, criteria);

    loss_struct_all{i} = loss_i;
    eff_struct_all{i} = eff_i;

    eff_scores(i) = min(struct2array(eff_i));  % conservative (min) score
end

%% 7. Sort and Display Top Designs
des_exact_rounding2(:, end+1) = eff_scores;
sorted_des = sortrows(des_exact_rounding2, size(des_exact_rounding2, 2), 'descend');

x_names = arrayfun(@(v) sprintf('%.0f', v), xx_temp, 'UniformOutput', false);
column_names = [x_names, {'Used Cost', 'Remaining', 'EffScore'}];
T = array2table(sorted_des, 'VariableNames', column_names);

disp('Top Rounded Designs (Multi-Objective under Budget C):');
disp(T(1:5, :));

%% 8. Merge Best Extended Design with Floor Allocation
y = [xx_temp; table2array(T(1, 1:end-3))];
all_keys = union(x(1,:), y(1,:));
combined_values = zeros(size(all_keys));

[~, loc_x] = ismember(x(1,:), all_keys);
[~, loc_y] = ismember(y(1,:), all_keys);

combined_values(loc_x) = combined_values(loc_x) + x(2,:);
combined_values(loc_y) = combined_values(loc_y) + y(2,:);

nz_idx = combined_values > 0;
design_filter = [all_keys(nz_idx); combined_values(nz_idx)];

%% 9. Final Evaluation: Efficiency
cxx_filter = 1 - q_cost + q_cost * design_filter(1,:);
w_final = design_filter(2,:) .* cxx_filter / C;

M_ex  = compute_FIM_GT_cost(design_filter(1,:), w_final, theta, q_cost);

eff_final = struct();
for i = 1:length(criteria)
    crit = criteria{i};
    switch crit
        case 'D'
            loss_ex = calc_loss_D(M_ex);
            eff_final.D = (loss_single.D / loss_ex)^(1/3);
        case 'A'
            loss_ex = calc_loss_A(M_ex);
            eff_final.A = loss_single.A / loss_ex;
        case 'c'
            loss_ex = calc_loss_c(M_ex, cVec_struct.cVec_c);
            eff_final.c = loss_single.c / loss_ex;
        case 'Ds'
            loss_ex = calc_loss_c(M_ex, cVec_struct.cVec_Ds);
            eff_final.Ds = loss_single.Ds / loss_ex;
        otherwise
            error('Unsupported criterion: %s', crit);
    end
end

%% 10. Final Output
disp("Before rounding");
disp(design_round1);

disp('Filtered Final Design (Exact Cost-Constrained):');
disp(design_filter);

fprintf('\nEfficiencies (Exact vs Approximate):\n');
criteria_list = fieldnames(eff_final);
for i = 1:length(criteria_list)
    crit = criteria_list{i};
    fprintf('  %s-efficiency: %.4f\n', crit, eff_final.(crit));
end

total_used_cost = sum(design_filter(2,:) .* (1 - q_cost + q_cost * design_filter(1,:)));
fprintf('\nTotal Budget: %.2f | Used: %.2f | Remaining: %.2f\n\n', ...
    C, total_used_cost, C - total_used_cost);