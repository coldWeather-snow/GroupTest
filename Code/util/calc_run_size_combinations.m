function results_sorted = calc_run_size_combinations(n_items, max_run_size)
    % n_items: number of item types
    % max_run_size: maximum allowed total number of items in a combination

    % Step 1: compute max counts per item type (since cost is 1, max is Cr)
    max_vals = repmat(max_run_size, 1, n_items);

    % Step 2: generate all combinations
    ranges = arrayfun(@(m) 0:m, max_vals, 'UniformOutput', false);
    [X{1:n_items}] = ndgrid(ranges{:});
    Xmat = cellfun(@(x) x(:), X, 'UniformOutput', false);
    combos = [Xmat{:}];

    % Step 3: compute total item count (run size)
    total_count = sum(combos, 2);

    % Step 4: filter combinations within the max run size
    valid_idx = total_count <= max_run_size;
    valid_combos = combos(valid_idx, :);
    valid_counts = total_count(valid_idx);
    remaining_slots = max_run_size - valid_counts;

    % Optional: filter to only full combinations (no remaining slots)
    % tight_idx = remaining_slots == 0;
    % results = [valid_combos(tight_idx, :), valid_counts(tight_idx)];

    results = [valid_combos, valid_counts, remaining_slots];

    % Step 5: sort by variable values (last to first)
    sort_columns = size(valid_combos, 2):-1:1;
    results_sorted = sortrows(results, sort_columns);
end