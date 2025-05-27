%% Compute OAD and OED for optimal group testing design using run-size-based rounding
clear;
start_time = tic;  % Start timer
addpath('../util');
criterion = 'Ds';
n = 9;               % total number of runs
M = 150;

%% Initialize parameters 
cVec_c = [0,1,1]';
cVec_Ds = [1,0,0]';
q_cost = 0;
S = [1, M]';
theta = [0.07, 0.93, 0.96]';
q = length(theta);
p0 = theta(1); p1 = theta(2); p2 = theta(3);
tol = 1E-5;
u = S(1):S(2);
N = length(u);
one_vec = ones(N,1); zero_vec = zeros(N,1);

%% CVX Optimization
cvx_begin quiet
cvx_precision best
variables w(N,1) del(1)
minimize del(1)
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
        disp("Error");
end
-w <= zero_vec;
one_vec' * w == 1;
cvx_end

kk = find(w > tol);
design_app = [u(kk); w(kk)'];
cvx_status

%% 3. Rounding Based on Run Size
xx = design_app(1,:);
ww = design_app(2,:);
n_prime = floor(n * ww);

deficit_indices = find(n_prime == 0);
excess_indices = find(n_prime >= 2);
for i = 1:length(deficit_indices)
    if isempty(excess_indices)
        error('Not enough room to adjust n_prime to make all entries >= 1.');
    end
    from_idx = excess_indices(1);
    to_idx = deficit_indices(i);
    n_prime(from_idx) = n_prime(from_idx) - 1;
    n_prime(to_idx) = n_prime(to_idx) + 1;
    if n_prime(from_idx) < 2
        excess_indices(1) = [];
    end
end
design_round1 = [xx; n_prime];
nr = round(n - sum(n_prime),4);

%% 4. Extended Support Points
n_index = 2;
result = [];
for i = 1:length(xx)
    result = [result, (xx(i)-n_index):(xx(i)+n_index)];
end
xx_temp = unique(result);
xx_temp = xx_temp(ismember(xx_temp, u));

%% 5. All Valid Run-Size Combinations
des_exact_rounding2 = calc_run_size_combinations(length(xx_temp), nr);

%% 6. Combine Approximate + Rounding Points
x = design_round1;
y = [xx_temp; des_exact_rounding2(:,1:end-2)];
all_keys = union(x(1,:), y(1,:));
n_cases = size(y,1)-1;
result_cell = cell(n_cases,1);

[~, loc_x] = ismember(x(1,:), all_keys);
[~, loc_y] = ismember(y(1,:), all_keys);

for k = 1:n_cases
    values = zeros(1,length(all_keys));
    values(loc_x) = values(loc_x) + x(2,:);
    values(loc_y) = values(loc_y) + y(k+1,:);
    result_cell{k} = [all_keys; values];
end

%% 7. Evaluate Designs
loss_round = zeros(n_cases,1);
for i = 1:n_cases
    des = result_cell{i};
    des(2,:) = des(2,:) / n;
    M = compute_FIM_GT_cost(des(1,:), des(2,:), theta, q_cost);
    switch criterion
        case 'D'
            loss_round_temp = calc_loss_D(M);
        case 'A'
            loss_round_temp = calc_loss_A(M);
        case 'c'
            loss_round_temp = calc_loss_c(M, cVec_c);
        case 'Ds'
            loss_round_temp = calc_loss_c(M, cVec_Ds);
        otherwise
            disp("Error");
    end
    loss_round(i) = loss_round_temp;
end

%% 8. Display Top Designs
des_exact_rounding2(:,end+1) = loss_round;
sorted_des = sortrows(des_exact_rounding2, size(des_exact_rounding2,2));

n_vars = size(sorted_des,2) - 3;
x_names = arrayfun(@(v) sprintf('%.0f', v), xx_temp, 'UniformOutput', false);
column_names = [x_names, {'Total Runs', 'Rem. Runs', 'Loss'}];
T = array2table(sorted_des, 'VariableNames', column_names);

disp('Original Design (Rounded):');
disp(design_round1);
disp('Top Rounded Designs (Run-size based):');
disp(T(1:5,:));

%% 9. Final Merging
x = design_round1;
y = [xx_temp; table2array(T(1,1:end-3))];
all_keys = union(x(1,:), y(1,:));
combined_values = zeros(size(all_keys));
[~, loc_x] = ismember(x(1,:), all_keys);
combined_values(loc_x) = combined_values(loc_x) + x(2,:);
[~, loc_y] = ismember(y(1,:), all_keys);
combined_values(loc_y) = combined_values(loc_y) + y(2,:);
nonzero_idx = combined_values > 0;
design_filter = [all_keys(nonzero_idx); combined_values(nonzero_idx)];

%% 10. Compute Efficiency
M_ex = compute_FIM_GT_cost(design_filter(1,:), design_filter(2,:)/n, theta, q_cost);
M_app = compute_FIM_GT_cost(design_app(1,:), design_app(2,:), theta, q_cost);

switch criterion
    case 'D'
        loss_ex = calc_loss_D(M_ex);
        loss_app = calc_loss_D(M_app);
        eff = (loss_app/loss_ex)^(1/3);
    case 'A'
        loss_ex = calc_loss_A(M_ex);
        loss_app = calc_loss_A(M_app);
        eff = loss_app / loss_ex;
    case 'c'
        loss_ex = calc_loss_c(M_ex, cVec_c);
        loss_app = calc_loss_c(M_app, cVec_c);
        eff = loss_app / loss_ex;
    case 'Ds'
        loss_ex = calc_loss_c(M_ex, cVec_Ds);
        loss_app = calc_loss_c(M_app, cVec_Ds);
        eff = loss_app / loss_ex;
    otherwise
        disp("Error");
end

%% 11. Output Results
disp('Final Combined Design:');
disp(design_filter);
disp('Approximate Design:');
disp(design_app);
fprintf('Efficiency: %.4f\n', eff);
fprintf('Total Runs: %d\n', n);
fprintf('Criterion: %s\n', criterion);

total_time = toc(start_time);  % End timer
fprintf('Total computational time: %.4f seconds\n', total_time);