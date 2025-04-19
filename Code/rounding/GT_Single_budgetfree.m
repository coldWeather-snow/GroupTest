% Description:
%   Compute OAD and OED for optimal group testing design using run-size-based rounding.

%% 1. Input
criterion = 'D';
calc_loss_D = @(my_FIM) 1/det(my_FIM);
calc_loss_D_logdet = @(my_FIM) -log(det(my_FIM));
calc_loss_A = @(my_FIM) trace(inv(my_FIM));
calc_loss_c = @(my_FIM, cc) cc' * inv(my_FIM) * cc;

cVec_c = [0,1,1]';
cVec_Ds = [1, 0,0]';
q_cost = 0;
S = [1, 150]'; % design space
theta = [0.07, 0.93, 0.96]';
q = length(theta);
p0 = theta(1); p1 = theta(2); p2 = theta(3);
n = 50;               % total number of runs (instead of budget)
tol = 1E-5;             % tolerance for filtering support points

%% 2. Initialization
u = S(1) : S(2);        % design space
N = length(u);
w = zeros(N+1, 1); del = 0;
one_vec = ones(N,1); zero_vec = zeros(N, 1);

%%% CVX Optimization
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

%%% Extract approximate optimal design
kk = find(w > tol);
design_app = [u(kk); w(kk)'];
cvx_status

%% 3. Rounding Based on Run Size
xx = design_app(1,:);
ww = design_app(2,:);
n_prime = floor(n * ww);

% Fix any zero entries by redistributing from indices with n ≥ 2
deficit_indices = find(n_prime == 0);
excess_indices = find(n_prime >= 2);

for i = 1:length(deficit_indices)
    if isempty(excess_indices)
        error('Not enough room to adjust n_prime to make all entries >= 1.');
    end

    from_idx = excess_indices(1);
    to_idx = deficit_indices(i);

    % Transfer one unit
    n_prime(from_idx) = n_prime(from_idx) - 1;
    n_prime(to_idx) = n_prime(to_idx) + 1;

    % If from_idx drops below 2, remove it from eligible donors
    if n_prime(from_idx) < 2
        excess_indices(1) = []; % remove current index
    end
end


design_round1 = [xx; n_prime];

nr = round(n - sum(n_prime), 4);  % remaining runs to allocate

%% 4. Generate Extended Support Points
n_index = 2;
result = [];

for i = 1:length(xx)
  startIdx = xx(i) - n_index;
  endIdx = xx(i) + n_index;
  result = [result, startIdx:endIdx];
end

xx_temp = unique(result);
xx_temp = xx_temp(ismember(xx_temp, u));

%% 5. Generate All Valid Run-Size Combinations
des_exact_rounding2 = calc_run_size_combinations(length(xx_temp), nr);  % NEW FUNCTION

%% 6. Combine Approximate + Rounding Design Points
x = design_round1;
y = [xx_temp; des_exact_rounding2(:, 1:end-2)];  % transpose combinations

% Unique keys from both x and y
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

%% 7. Evaluate Designs via Fisher Information Matrix
loss_round = zeros(n_cases, 1);
for i = 1:n_cases
  des = result_cell{i};
  des(2, :) = des(2, :) / n;  % normalize by total runs
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
des_exact_rounding2(:, end+1) = loss_round;
sorted_des = sortrows(des_exact_rounding2, size(des_exact_rounding2, 2));

n_vars = size(sorted_des, 2) - 3;
x_names = arrayfun(@(v) sprintf('%.0f', v), xx_temp, 'UniformOutput', false);
column_names = [x_names, {'Total Runs', 'Rem. Runs', 'Loss'}];

T = array2table(sorted_des, 'VariableNames', column_names);

disp('Original Design (Rounded):');
disp(design_round1);

disp('Top Rounded Designs (Run-size based):');
disp(T(1:5, :));

%%% Algin the output to print better
% Display Top 5 Rounded Designs Nicely
fprintf('\nTop Rounded Designs (Run-size based):\n\n');

% Extract top 5
T_top = T(1:5, :);


design_round1

%%
x = design_round1;
y = [xx_temp; table2array(T(1, 1:end-3))];

% Step 1: Get unique keys
all_keys = union(x(1,:), y(1,:));

% Step 2: Initialize result vector
combined_values = zeros(size(all_keys));

% Step 3: Map and add values from x
[~, loc_x] = ismember(x(1,:), all_keys);
combined_values(loc_x) = combined_values(loc_x) + x(2,:);

% Step 4: Map and add values from y
[~, loc_y] = ismember(y(1,:), all_keys);
combined_values(loc_y) = combined_values(loc_y) + y(2,:);

% Step 5: Filter out zero weights
nonzero_idx = combined_values > 0;
filtered_keys = all_keys(nonzero_idx);
filtered_values = combined_values(nonzero_idx);

% Step 6: Return final design
design_filter = [filtered_keys; filtered_values];


M_ex = compute_FIM_GT_cost(design_filter(1,:), design_filter(2,:)/n, theta, q_cost);
M_app = compute_FIM_GT_cost(design_app(1,:), design_app(2,:), theta, q_cost);
 switch criterion
    case 'D'
      loss_ex = calc_loss_D(M_ex);
      loss_app = calc_loss_D(M_app);
      eff = (loss_app / loss_ex )^(1/3);
    case 'A'
      loss_ex = calc_loss_A(M_ex);
      loss_app = calc_loss_A(M_app);
      eff =  loss_app / loss_ex;
    case 'c'
      loss_ex = calc_loss_c(M_ex, cVec_c);
      loss_app = calc_loss_c(M_app, cVec_c);
      eff =  loss_app / loss_ex;
    case 'Ds'
      loss_ex = calc_loss_c(M_ex, cVec_Ds);
      loss_app = calc_loss_c(M_app, cVec_Ds);
      eff =  loss_app / loss_ex;
    otherwise
      disp("Error");
 end

 design_filter
 design_app

eff
n
criterion