function eff_struct = compute_efficiencies(loss_single, loss_target, criteria)
% COMPUTE_EFFICIENCIES Computes efficiencies for multiple criteria.
%
% Inputs:
%   loss_single     - struct with reference (optimal) losses (e.g., from D-opt design)
%   loss_target  - struct with losses of current design (e.g., from maximin)
%   criteria     - cell array of strings, e.g., {'D', 'A', 'Ds'}
%
% Output:
%   eff_struct   - struct with efficiency values for each criterion

    eff_struct = struct();

    for i = 1:numel(criteria)
        crit = criteria{i};

        if ~isfield(loss_single, crit) || ~isfield(loss_target, crit)
            error('Missing loss value for criterion: %s', crit);
        end

        switch crit
            case 'D'
                % D-efficiency is typically: (|M_opt| / |M_target|)^(1/q)
                % But if losses are 1/|M|, then:
                eff_struct.D = (loss_single.D / loss_target.D)^(1/3);  % 3 = dim(FIM)
            otherwise
                % For A, Ds, c: efficiency = L_opt / L_target
                eff_struct.(crit) = loss_single.(crit) / loss_target.(crit);
        end
    end
end