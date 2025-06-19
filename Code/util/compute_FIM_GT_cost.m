function [FIM, X] = compute_FIM_GT_cost(di, wi, theta, q_cost)
    % COMPUTE_FIM - Computes the Fisher Information Matrix (FIM) in a vectorized manner.
    %
    % Inputs:
    %   di   - Column vector of design support points.
    %   wi   - Column vector of corresponding weights.
    %   theta - Probability parameters [p0, p1, p2].
    %
    % Outputs:
    %   FIM  - Computed Fisher Information Matrix (3 × 3).
    
    % Extract probability parameters
    p0 = theta(1);
    p1 = theta(2);
    p2 = theta(3);

    % Ensure column vector format
    x_vals = di(:);
    w_vals = wi(:);
    c_vals = (1 - q_cost) + q_cost .* x_vals ;

    % Compute probability values (Vectorized)
    pi_vals = p1 - (p1 + p2 - 1) * (1 - p0) .^ x_vals;
    lam_vals = 1 ./ (c_vals .* pi_vals .* (1 - pi_vals));  % Element-wise inverse

    % Compute f vectors for all x values (Vectorized)
    f0_vals = x_vals .* (p1 + p2 - 1) .* (1 - p0) .^ (x_vals - 1);
    f1_vals = 1 - (1 - p0) .^ x_vals;
    f2_vals = -(1 - p0) .^ x_vals;

    % Stack f vectors as a matrix (3 × N)
    F_vals = [f0_vals'; f1_vals'; f2_vals'];  % Transpose to match dimensions
    X = F_vals';
    % Compute weighted Fisher Information Matrix (Vectorized)
    FIM = F_vals * diag(lam_vals .* w_vals) * F_vals';
end