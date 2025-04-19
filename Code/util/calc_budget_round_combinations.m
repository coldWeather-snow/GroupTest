function results_sorted = calc_budget_round_combinations(cxx, Cr)
    n = length(cxx);  % number of item types

    % Step 1: compute max counts
    max_vals = floor(Cr ./ cxx);

    % Step 2: generate all combinations
    ranges = arrayfun(@(m) 0:m, max_vals, 'UniformOutput', false);
    [X{1:n}] = ndgrid(ranges{:});
    Xmat = cellfun(@(x) x(:), X, 'UniformOutput', false);
    combos = [Xmat{:}];

    % Step 3: compute total cost and remainder with rounding to prevent float overflow
    total_cost = round(combos * cxx(:), 6);  % round to 6 decimals
    Cr_rounded = round(Cr, 6);               % also round Cr for safe comparison
    remaining_budget = round(Cr_rounded - total_cost, 6);

    valid_idx = total_cost <= Cr_rounded;
    valid_combos = combos(valid_idx, :);
    valid_costs = total_cost(valid_idx);
    remaining_budget = remaining_budget(valid_idx);

    % Step 4: filter for minimal remaining budget (less than min item cost)
    min_item_cost = min(cxx);
    tight_idx = remaining_budget < min_item_cost - 1e-6;  % small tolerance
    results = [valid_combos(tight_idx, :), valid_costs(tight_idx), remaining_budget(tight_idx)];

    % Step 5: sort by variable values (last to first)
    sort_columns = size(valid_combos, 2):-1:1;
    results_sorted = sortrows(results, sort_columns);
end