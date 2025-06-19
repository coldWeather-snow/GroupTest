function eta = calc_eta_weights(tstar, loss_ref, loss_model, directional_derivatives, criterion, tol)
% SOLVE_ETA_WEIGHTS solves for optimal eta weights for multi-objective optimal design.
%
% Inputs:
%   tstar                   - scalar (maximin t*)
%   loss_ref               - struct with fields like 'D', 'A', 'Ds', 'c'
%   loss_model             - struct with same fields (actual model loss values)
%   directional_derivatives - struct with fields: 'dD', 'dA', etc., each Nx1
%   criterion              - cell array of strings, e.g., {'D', 'A', 'Ds'}
%   tol                    - scalar tolerance for constraint (e.g., 1e-5)
%
% Output:
%   eta                    - vector of optimal weights (K x 1)

    K = numel(criterion);
    N = length(directional_derivatives.(sprintf('d%s', criterion{1})));  % assume all same length

    cvx_begin quiet
        cvx_precision best
        variable eta(K, 1)
        minimize sum(eta)

        % === Constraint 1: Normalization ===
        norm_expr = 0;
        for k = 1:K
            name = criterion{k};
            if strcmp(name, 'D')
                norm_expr = norm_expr + eta(k) * (3 / tstar);  % D-opt uses log det
            elseif strcmp(name, 'E')
                norm_expr = norm_expr + eta(k) * -loss_ref.E/tstar^2;  % D-opt uses log det
            else
                norm_expr = norm_expr + eta(k) * loss_ref.(name);  % others use linear loss
            end
        end
        norm_expr == 1;

        % === Constraint 2: Loss closeness constraints ===
        for k = 1:K
            name = criterion{k};
            switch name
                case 'D'
                    delta = log(loss_model.D) - log(loss_ref.D) - 3 * log(tstar);
              case 'E'
                  delta = loss_model.E - loss_ref.E/tstar;
                otherwise
                    delta = loss_model.(name) - tstar * loss_ref.(name);
            end
            eta(k) * delta <= tol;
            -eta(k) * delta <= tol;
        end

        % === Constraint 3: Directional derivatives pointwise ≤ tol ===
        sum_deriv = zeros(N, 1);
        for k = 1:K
            field = ['d' criterion{k}];
            sum_deriv = sum_deriv + eta(k) * directional_derivatives.(field);
        end
        sum_deriv <= tol;

        % === Constraint 4: Non-negativity ===
        -eta <= 0;
    cvx_end
end