function eff_struct = compute_efficiencies(loss_single, loss_target, criteria)
% COMPUTE_EFFICIENCIES Computes efficiencies for multiple criteria.
%
% Inputs:
%   loss_single - struct with optimal losses (e.g., from D-opt design)
%   loss_target - struct with losses of current design
%   criteria    - string or cell array of strings, e.g., 'D' or {'D','A'}
%
% Output:
%   eff_struct  - struct with efficiency values

    if ischar(criteria) || isstring(criteria)
        criteria = {char(criteria)};  % wrap in cell array
    end

    eff_struct = struct();

    for i = 1:numel(criteria)
        crit = criteria{i};

        if ~isfield(loss_single, crit) || ~isfield(loss_target, crit)
            error('Missing loss value for criterion: %s', crit);
        end

        switch crit
            case 'D'
                eff_struct.D = (loss_single.D  / loss_target.D )^(1/3);  % 3 = dim(FIM)
            case 'E'
                eff_struct.E = loss_target.E / loss_single.E;
            otherwise
                eff_struct.(crit) = loss_single.(crit) / loss_target.(crit);
        end
    end
end