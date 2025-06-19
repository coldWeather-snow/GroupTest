function [lambda_min_E, r_star, Q] = calc_smallest_eig_info(M, tol)
%GET_SMALLEST_EIG_INFO Extracts smallest eigenvalue info for symmetric matrix
%
% Inputs:
%   M   : Symmetric positive semidefinite matrix (q × q)
%   tol : Tolerance for eigenvalue multiplicity check (default = 1e-8)
%
% Outputs:
%   lambda_min_E : The smallest eigenvalue of M
%   r_star       : Multiplicity of the smallest eigenvalue
%   Q            : q × r_star matrix of orthonormal eigenvectors corresponding
%                  to the smallest eigenvalue

    if nargin < 2
        tol = 1e-8;
    end

    [V, D] = eig(M);
    eigenvalues = diag(D);
    lambda_min_E = min(eigenvalues);

    % Multiplicity and eigenvectors
    r_star = sum(abs(eigenvalues - lambda_min_E) < tol);
    Q = V(:, abs(eigenvalues - lambda_min_E) < tol);
end