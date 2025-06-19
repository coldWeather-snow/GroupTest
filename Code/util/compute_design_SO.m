function result = compute_design_SO(u, theta, q_cost, q, criterion, opts, tol)
% COMPUTE_OPTIMAL_DESIGN
% Computes the optimal design for a given criterion in group testing.
%
% Inputs:
%   u         - vector of candidate design points
%   theta     - model parameters [p0; p1; p2]
%   q_cost    - cost parameter
%   q         - dimension of the FIM (should match length(theta))
%   criterion - string: 'D', 'A', 'Ds', 'E' or 'c'
%   opts      - struct with optional vectors:
%                 .cVec   (for 'c')
%                 .DsVec  (for 'Ds')
%   tol       - tolerance to filter zero weights
%
% Output:
%   result    - struct with fields:
%                 .criterion
%                 .design (2 x k)
%                 .loss
%                 .w (N x 1)
%                 .M (FIM matrix)

N = length(u);
one_vec = ones(N, 1);
zero_vec = zeros(N, 1);

cvx_begin quiet
cvx_precision best
variables w(N,1) del(1)
expression M(q, q)
minimize(del)
M = compute_FIM_GT_cost(u, w, theta, q_cost);
switch criterion
  case 'D'
    -log_det(M) <= del
  case 'A'
    trace_inv(M) <= del;
  case 'Ds'
    if ~isfield(opts, 'cVec_Ds')
      error('Missing DsVec in opts for Ds-optimality.');
    end
    cVec = opts.cVec_Ds;

    subject to
    matrix_frac(cVec, M) <= del;
  
  case 'c'
    if ~isfield(opts, 'cVec_c')
      error('Missing cVec in opts for c-optimality.');
    end
    cVec = opts.cVec_c;
    matrix_frac(cVec, M) <= del;

  case 'E'
    -lambda_min(M) < del;
  otherwise
    error('Unsupported criterion: %s', criterion);
end
subject to
-w <= zero_vec;
one_vec' * w == 1;
cvx_end

% Filter support points
idx = find(w > tol);
design = [u(idx); w(idx)'];

% Compute loss
switch criterion
  case 'D'
    loss = calc_loss_D(M);
  case 'A'
    loss = calc_loss_A(M);
  case {'c', 'Ds'}
    loss = calc_loss_c(M, cVec);
  case 'E'
    loss = calc_loss_E(M);
end

% Return all outputs in a struct
result = struct( ...
  'M', M, ...
  'criterion', criterion, ...
  'design', design, ...
  'loss', loss, ...
  'cvx_status', cvx_status...
  );
end