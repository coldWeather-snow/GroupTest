function result = compute_maximin_design(u, theta, q_cost, q, loss_ref, criteria, cVec_struct, tol)
% COMPUTE_MAXIMIN_DESIGN solves a Maximin optimal design problem.
%
% Inputs:
%   u         - vector of design points
%   theta     - parameter vector
%   q_cost    - cost coefficient
%   q         - dimension of FIM
%   loss_ref  - struct with reference loss values (e.g., from single-optimality designs)
%   criteria  - cell array of strings, e.g., {'D', 'A', 'Ds', 'c'}
%   cVec_struct      - struct with .cVec, .DsVec if needed
%   tol       - tolerance for weight filtering
%
% Output:
%   result    - struct with fields:
%                .criterion
%                .design
%                .loss
%                .efficiency
%                .w
%                .tstar
%                .M

    N = length(u);
    % Build CVX optimization
    cvx_begin quiet
    cvx_precision best
    variables w(N,1) tstar(1)
    expression M(q, q)
    M = compute_FIM_GT_cost(u, w, theta, q_cost);
    minimize tstar

    % General constraints
    -tstar <= 0;
    0 <= w;
    sum(w) == 1;

    % Criterion-specific constraints
    for k = 1:length(criteria)
        crit = criteria{k};
        switch crit
            case 'D'
                -log_det(M) - (log(loss_ref.D) + q * log(tstar)) <= 0;
            case 'A'
                trace_inv(M) - tstar * loss_ref.A <= 0  ;
            case 'Ds'
                if ~isfield(cVec_struct, 'cVec_Ds')
                    error('Missing DsVec in opts for Ds-optimality.');
                end
                matrix_frac(cVec_struct.cVec_Ds, M) - tstar * loss_ref.Ds <= 0;
            case 'c'
                if ~isfield(cVec_struct, 'cVec_c')
                    error('Missing cVec in opts for c-optimality.');
                end
                matrix_frac(cVec_struct.cVec_c, M) - tstar * loss_ref.c <= 0 ; 
          case 'E'
            -lambda_min(M) - loss_ref.E * inv_pos(tstar) <= 0 ;
            otherwise
                error('Unsupported criterion: %s', crit);
        end
    end
cvx_end

    % Extract support design
    idx = find(w > tol);
    design = [u(idx); w(idx)'];

    % Compute actual losses and efficiencies
    loss = compute_losses(M, criteria, cVec_struct);
    eff = compute_efficiencies(loss_ref, loss, criteria);

    result = struct( ...
        'criterion', {criteria}, ...
        'design', design, ...
        'loss', loss, ...
        'efficiency', eff, ...
        'tstar', tstar, ...
        'M', M...
    );
end