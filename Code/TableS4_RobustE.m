%Compute the maximin design with dual objective for group testing
%experiment

%% addpath
% clear;
start_time = tic;  % Start timer
addpath('util');

rho = 0.01 % uncertainty

%% Set up the parameters
q_cost = 0.8;
% tol = 1E-5;
M = 61;
S = [1, M]'; % design space
cVec_c = [0, 1, 1]';
cVec_Ds = [1, 0, 0]';
cVec_struct = struct('cVec_c', cVec_c, 'cVec_Ds', cVec_Ds);
theta = [0.07, 0.93, 0.96]';
p0 = theta(1); p1 = theta(2); p2 = theta(3);
q = length(theta);
u =  S(1) : S(2);
N = length(u);

% q = p

cvx_begin sdp; cvx_precision best    
    variable w(N) nonnegative 
    variables t(1) s(1) 

    [FIM,V] = compute_FIM_GT_cost(u, w, theta, q_cost);
    minimize -s
    subject to
            FIM - rho * sqrt(N) *t* eye(q) >= s * eye(q)
        norm(w, 2) <= t
        w >= 0; sum(w) == 1;
cvx_end

loss_rho = cvx_optval;


% --- Step 4: Extract approximate support points ---
tol = 1e-4;
idx = find(w > tol);
design = [u(idx); w(idx)']
round(design, 3)


%% Robust design at rho = 0, the regular E-opt
rho_0 = 0;
cvx_begin sdp; cvx_precision best    
    variable w(N) nonnegative 
    variables t_rho0(1) s_rho0(1) 

    [FIM_rho0,V_rho0] = compute_FIM_GT_cost(u, w, theta, q_cost);
    minimize -s_rho0
    subject to
            FIM_rho0 - rho_0 * sqrt(N) *t_rho0* eye(q) >= s_rho0 * eye(q)
        norm(w, 2) <= t_rho0
        w >= 0; sum(w) == 1;
cvx_end

loss_rho0 = cvx_optval;
% --- Step 4: Extract approximate support points ---
tol = 1e-4;
idx = find(w > tol);
design_rho0 = [u(idx); w(idx)']
round(design_rho0, 3)

design = [u(idx); w(idx)']
round(design, 3)

lambda_min(FIM)
lambda_min(FIM_rho0)
-s
-s_rho0

round(loss_rho/loss_rho0,3)
