function d_out = calc_directional_derivatives(u, M, theta, q_cost, criterion)
% Compute directional derivatives for given criteria
%
% Inputs:
%   u         - design space (vector of integers)
%   M         - Fisher Information Matrix (FIM)
%   p0, p1, p2- model parameters
%   q_cost    - cost coefficient
%   criterion - cell array of strings: e.g., {'D', 'A', 'Ds'}
%
% Output:
%   d_out - struct with fields dD, dA, dDs depending on criterion
p0  = theta(1);
p1 = theta(2);
p2 = theta(3);
N = length(u);
q = size(M,1);
inv_M = inv(M);
tr_inv_M = trace(inv_M);

% Preallocate outputs if requested
d_out = struct();
if any(strcmp(criterion, 'D'))
  d_out.dD = zeros(N, 1);
end
if any(strcmp(criterion, 'A'))
  d_out.dA = zeros(N, 1);
end
if any(strcmp(criterion, 'Ds'))
  d_out.dDs = zeros(N, 1);
  cVec_Ds = [1; 0; 0]; % This can be passed as an optional param if needed
  cVec_Ds_Minv = inv_M * cVec_Ds;
  cVec_Ds_norm = cVec_Ds' * inv_M * cVec_Ds;
end

for i = 1:N
  x = u(i);
  % Compute f(x)
  fx0 = x*(p1 + p2 -1) * (1-p0)^(x-1);
  fx1 = 1 - (1-p0)^x;
  fx2 = - (1 - p0)^x;
  fx = [fx0; fx1; fx2];

  % Compute λ(x)
  cx = 1 - q_cost + q_cost * x;
  pix = p1 - (p1 + p2 -1) * (1-p0)^x;
  lambda_x = 1 / (cx * pix * (1 - pix));

  fx_fxT = fx * fx';

  % Compute requested directional derivatives
  if isfield(d_out, 'dD')
    d_out.dD(i) = lambda_x * trace(inv_M * fx_fxT) - q;
  end
  if isfield(d_out, 'dA')
    d_out.dA(i) = lambda_x * trace(inv_M * fx_fxT * inv_M) - tr_inv_M;
  end
  if isfield(d_out, 'dDs')
    d_out.dDs(i) = lambda_x * (cVec_Ds' * inv_M * fx_fxT * inv_M * cVec_Ds) - cVec_Ds_norm;
  end
end
end