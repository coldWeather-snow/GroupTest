%% Run-Size-Based Rounding for Multi-Objective Design
% This script rounds an approximate optimal design to exact run counts.
% Assumes q_cost = 0 and one criterion is used to evaluate efficiency.

%% 1. Setup
n = 100;                        % Total number of runs

% Load or define design_struct (e.g., from compute_maximin_design)


design_app = result_DADs;      % input
criteria = {'D', 'A', 'Ds'}



% design_app = result_DDsc;      % input
% criteria = {'D',  'Ds', 'c'};  

%% 2. First Round: Integer Approximation
n_prime = floor(n * design_app.design(2,:));

% Fix 0 counts by redistributing from n ≥ 2 entries
deficit_idx = find(n_prime == 0);
excess_idx  = find(n_prime >= 2);

for i = 1:length(deficit_idx)
    if isempty(excess_idx)
        error('Not enough excess runs to fix zero entries.');
    end
    from = excess_idx(1);
    to = deficit_idx(i);

    n_prime(from) = n_prime(from) - 1;
    n_prime(to) = n_prime(to) + 1;

    if n_prime(from) < 2
        excess_idx(1) = [];
    end
end

design_round1 = [design_app.design(1,:); n_prime];
nr = round(n - sum(n_prime), 4);  % leftover runs

%% 3. Extended Support Points
n_index = 2;
extended_pts = [];
xx = design_round1(1,:);
ww = design_round1(2,:);

for i = 1:length(xx)
    extended_pts = [extended_pts, (xx(i) - n_index):(xx(i) + n_index)];
end

xx_temp = unique(extended_pts);
xx_temp = xx_temp(ismember(xx_temp, u));

%% 4. Generate Valid Run Combinations
des_exact_rounding2 = calc_run_size_combinations(length(xx_temp), nr);  % External function

%% 5. Combine Rounded + Extended Combinations
x = design_round1;
y = [xx_temp; des_exact_rounding2(:, 1:end-2)];

all_keys = union(x(1,:), y(1,:));
n_cases = size(y, 1) - 1;
result_cell = cell(n_cases, 1);

[~, loc_x] = ismember(x(1,:), all_keys);
[~, loc_y] = ismember(y(1,:), all_keys);

for k = 1:n_cases
    values = zeros(1, length(all_keys));
    values(loc_x) = values(loc_x) + x(2,:);
    values(loc_y) = values(loc_y) + y(k+1, :);
    result_cell{k} = [all_keys; values];
end

%% 6. Evaluate Designs: multi-criteria loss + efficiency


loss_struct_all = cell(n_cases, 1);
eff_struct_all = cell(n_cases, 1);
eff_scores = zeros(n_cases, 1);  % for ranking

for i = 1:n_cases
    des = result_cell{i};
    des = des(:, des(2,:) > 0);
    nnnn = sum(des(2,:));

    % des(2,:) = des(2,:) / nnnn;  % Normalize to proportions
    des(2,:) = des(2,:) / n;  % Normalize to proportions
    
    M = compute_FIM_GT_cost(des(1,:), des(2,:), theta, q_cost);

    % Compute multi-criterion loss
    loss_i = compute_losses(M, criteria, cVec_struct);
    
    % Compute efficiencies vs single-objective benchmarks
    eff_i = compute_efficiencies(loss_single, loss_i, criteria);
    % Store
    loss_struct_all{i} = loss_i;
    eff_struct_all{i} = eff_i;

    % Score for sorting (mean or geomean works)
    eff_scores(i) = min(struct2array(eff_i));
end

%% 7. Sort and Tabulate Top Designs
% Append efficiency score to run combination table
des_exact_rounding2(:, end+1) = eff_scores;

% Sort by efficiency score (descending)
sorted_des = sortrows(des_exact_rounding2, size(des_exact_rounding2, 2), 'descend');

% Column names
x_names = arrayfun(@(v) sprintf('%.0f', v), xx_temp, 'UniformOutput', false);
column_names = [x_names, {'Total Runs', 'Rem. Runs', 'AvgEff'}];

% Create table
T = array2table(sorted_des, 'VariableNames', column_names);

% Display top 5 designs
disp('Top Rounded Designs:');
disp(T(1:5, :));

%% 8. Merge Best Extension with Rounded Design
y = [xx_temp; table2array(T(1, 1:end-3))];
all_keys = union(x(1,:), y(1,:));
combined_values = zeros(size(all_keys));

[~, loc_x] = ismember(x(1,:), all_keys);
[~, loc_y] = ismember(y(1,:), all_keys);

combined_values(loc_x) = combined_values(loc_x) + x(2,:);
combined_values(loc_y) = combined_values(loc_y) + y(2,:);

nz_idx = combined_values > 0;
design_filter = [all_keys(nz_idx); combined_values(nz_idx)];

%% 9. Compute Final Efficiency for All Criteria

% Compute Fisher information matrices
M_ex = compute_FIM_GT_cost(design_filter(1,:), design_filter(2,:) / n, theta, q_cost);
M_app = compute_FIM_GT_cost(design_app.design(1,:), design_app.design(2,:), theta, q_cost);

eff_final = struct();  % store all efficiencies

for i = 1:length(criteria)
    crit = criteria{i};

    switch crit
        case 'D'
            loss_ex = calc_loss_D(M_ex);
            
            eff_final.D = (loss_single.D/ loss_ex)^(1/3);  % corrected

        case 'A'
            loss_ex = calc_loss_A(M_ex);
            
            eff_final.A = loss_single.A / loss_ex;  % corrected

        case 'c'
            loss_ex = calc_loss_c(M_ex, cVec_struct.cVec_c);
            
            eff_final.c = loss_single.c / loss_ex;  % corrected

        case 'Ds'
            loss_ex = calc_loss_c(M_ex, cVec_struct.cVec_Ds);
            
            eff_final.Ds = loss_single.Ds/ loss_ex;  % corrected

        otherwise
            error('Unsupported criterion: %s', crit);
    end
end
%% 10. Display Final Result
disp('Filtered Final Design (Exact Run-Sized):');
disp(design_filter);

fprintf('\nEfficiencies (Exact vs Approximate) by Criterion:\n');
criteria = fieldnames(eff_final);

for i = 1:length(criteria)
    crit = criteria{i};
    fprintf('  %s-efficiency: %.4f\n', crit, eff_final.(crit));
end

fprintf('\nTotal Runs Allocated: %d\n', n);

nr

design_round1