function loss_struct = compute_losses(M, criteria, opts)
% COMPUTE_LOSSES Computes losses for specified optimality criteria
%
% Inputs:
%   M         - Fisher Information Matrix (q x q)
%   criteria  - Cell array of criterion strings, e.g., {'D', 'A', 'c'}
%   opts      - (Optional) Struct with fields: .cVec_c, .cVec_Ds
%
% Output:
%   loss_struct - struct with fields like .D, .A, .c, .Ds

    if nargin < 3
        opts = struct();  % Use empty struct if not provided
    end

    loss_struct = struct();

    for i = 1:numel(criteria)
        crit = criteria{i};
        switch crit
            case 'D'
                loss_struct.D = calc_loss_D(M);

            case 'A'
                loss_struct.A = calc_loss_A(M);

            case 'c'
                if ~isfield(opts, 'cVec_c')
                    error('Missing cVec in opts for c-optimality.');
                end
                loss_struct.c = calc_loss_c(M, opts.cVec_c);

            case 'Ds'
                if ~isfield(opts, 'cVec_Ds')
                    error('Missing cVec_Ds in opts for Ds-optimality.');
                end
                loss_struct.Ds = calc_loss_c(M, opts.cVec_Ds);

            otherwise
                error('Unsupported criterion: %s', crit);
        end
    end
end