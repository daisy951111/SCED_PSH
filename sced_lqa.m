% File: /Users/zliu20/Documents/Playground/hlselo_simulation_improved.m
function out = sced_lqa(cfg)
% Improved H-LSELO Monte Carlo simulation
% - Fixes covariance and indexing issues
% - Vectorizes heavy loops
% - Adds reproducible experiment wrapper and file exports

if nargin < 1
    cfg = default_cfg();
else
    cfg = apply_defaults(cfg, default_cfg());
end

rng(cfg.seed, 'twister');
tic;

% Core parameters
n = cfg.n;
p = cfg.p;
N = cfg.N;
mu = zeros(p, 1);
a = 1:p;
sigma = cfg.rho .^ abs(a - a');

tk_ini = cfg.tkini;
tk_main = cfg.tkmain;

% % True beta
% if strcmpi(cfg.signal_case, 'strong')
%     Beta = [-1.0; 1.0; 0; 0; 0; -1.0; zeros(p - 6, 1)];
% else
%     Beta = 0.5 * [-1.0; 1.0; 0; 0; 0; -1.0; zeros(p - 6, 1)];
% end
% index = find(Beta ~= 0);

% True beta: user-specified cfg.Beta takes priority
if isfield(cfg, 'Beta') && ~isempty(cfg.Beta)
    Beta = cfg.Beta(:);
    if numel(Beta) ~= p
        error('cfg.Beta must have length p.');
    end
else
    if strcmpi(cfg.signal_case, 'strong')
        Beta = [-1.0; 1.0; 0; 0; 0; -1.0; zeros(p - 6, 1)];
    else
        Beta = 0.5 * [-1.0; 1.0; 0; 0; 0; -1.0; zeros(p - 6, 1)];
    end
end
index = find(Beta ~= 0);

% Tuning grids
if isempty(cfg.theta)
    if strcmpi(cfg.signal_case, 'strong')
        theta = (0.60:-0.005:0.50) * 1.25e4;
    else
        theta = (0.575:-0.002:0.535) * 1.25e4;
    end
else
    theta = cfg.theta;
end
if isempty(cfg.delta)
    if strcmpi(cfg.signal_case, 'strong')
        delta = 1.5;
    else
        delta = 0.057;
    end
else
    delta = cfg.delta;
end

% Storage
Censorrate = zeros(N, 1);
initial_beta = zeros(p, N);
hlselo = zeros(p, N);
opt_theta = zeros(1, N);
opt_delta = zeros(1, N);
se_hlselo = zeros(N, 1);
cov_beta = zeros(numel(index), N);

Q = cfg.Q;
initial_beta_R = zeros(p, Q);

for iter = 1:N
    fprintf('iter = %d\n', iter);
    rng(iter, 'twister');

    % robust initial values via averaging Q replicates
    for j = 1:Q
        seed_ij = j * N * iter;
        rng(seed_ij, 'twister');
        [Z0, X0, T0, C0, Iota0, Delta0, R0] = survival_data(n, Beta, mu, sigma, cfg.tau, seed_ij);
        W0 = Weight(X0, T0, C0, Delta0, n);
        initial_beta_R(:, j) = ini_beta(Z0, cfg.tol, Delta0, Iota0, R0, W0, tk_ini);
    end
    initial_beta(:, iter) = mean(initial_beta_R, 2);

    % simulation dataset
    [Z, X, T, C, Iota, Delta, R] = survival_data(n, Beta, mu, sigma, cfg.tau, iter);
    Censorrate(iter) = 1 - mean(Delta);
    W = Weight(X, T, C, Delta, n);

    [hlselo(:, iter), opt_theta(iter), opt_delta(iter)] = hlselo_lqa(...
        n, initial_beta(:, iter), Z, cfg.tol, Delta, Iota, R, W, theta / 0.5, delta, tk_main);

    se_hlselo(iter) = (hlselo(:, iter) - Beta)' * sigma * (hlselo(:, iter) - Beta);
    cov_beta(:, iter) = cov_PSH(index, Delta, Iota, Z, W, R, opt_delta(iter), opt_theta(iter), hlselo(:, iter));
end

% Assessment
BIAS = mean(hlselo(index, :), 2) - Beta(index);
SE = std(hlselo(index, :), 0, 2);
ASE = mean(sqrt(max(cov_beta, 0)), 2);
ASE_std = std(sqrt(max(cov_beta, 0)), 0, 2);

lower = hlselo(index, :) - 1.96 * sqrt(max(cov_beta, 0));
upper = hlselo(index, :) + 1.96 * sqrt(max(cov_beta, 0));
in_interval = (lower <= Beta(index)) & (upper >= Beta(index));
CP_beta = mean(in_interval, 2);

% Exact model-selection correctness
active_true = all(hlselo(index, :) ~= 0, 1);
inactive_false = all(hlselo(setdiff(1:p, index), :) == 0, 1);
corr = mean(active_true & inactive_false);

MSE = mean(se_hlselo);
N_plus = sum(sum(hlselo(setdiff(1:p, index), :) ~= 0)) / N;
N_minus = sum(sum(hlselo(index, :) == 0)) / N;
Size = sum(sum(hlselo ~= 0)) / N;

% Tables
param_names = compose('beta_%d', index);
T = table(param_names(:), BIAS, SE, ASE, ASE_std, CP_beta, ...
    'VariableNames', {'Parameters', 'BIAS', 'SE', 'ASE', 'ASE_std', 'CP_beta'});

T1 = table("HSELO", corr, MSE, N_plus, N_minus, Size, ...
    'VariableNames', {'Method', 'Pcorr', 'MSE', 'FPR', 'PNR', 'Size'});

CR = mean(Censorrate);
run_time = toc;

out = struct();
out.config = cfg;
out.Beta = Beta;
out.index = index;
out.Censorrate = Censorrate;
out.initial_beta = initial_beta;
out.hlselo = hlselo;
out.opt_theta = opt_theta;
out.opt_delta = opt_delta;
out.cov_beta = cov_beta;
out.metrics = struct('CR', CR, 'runtime_sec', run_time, 'Pcorr', corr, 'MSE', MSE, ...
    'FPR', N_plus, 'PNR', N_minus, 'Size', Size);
out.table_coef = T;
out.table_model = T1;

if cfg.verbose
    disp(T);
    disp(T1);
    fprintf('CR = %.4f\n', CR);
    fprintf('time = %.2f sec\n', run_time);
end

if cfg.save_results
    if ~exist(cfg.out_dir, 'dir')
        mkdir(cfg.out_dir);
    end

    tag = sprintf('hlselo_n%d_p%d_N%d_rho%.2f_%s_seed%d', ...
    n, p, N, cfg.rho, builtin('lower', cfg.signal_case), cfg.seed);

    % tag = sprintf('hlselo_n%d_p%d_N%d_rho%.2f_%s_seed%d', ...
    %     n, p, N, cfg.rho, lower(cfg.signal_case), cfg.seed);

    save(fullfile(cfg.out_dir, [tag '.mat']), 'out');
    writetable(T, fullfile(cfg.out_dir, [tag '_coef.csv']));
    writetable(T1, fullfile(cfg.out_dir, [tag '_model.csv']));
end
end


function cfg = default_cfg()
cfg = struct();
cfg.n = 300;
cfg.p = 10;
cfg.N = 200;
cfg.rho = 0.25;
cfg.rho_set = [0.25, 0.5];
cfg.signal_case = 'weak'; % 'weak' or 'strong'
cfg.signal_cases = {'weak', 'strong'};
cfg.auto_tau = true;
cfg.tau = 1.60;
cfg.theta = [];
cfg.delta = [];
cfg.Q = 2;
cfg.tkini = 4;   % step size in ini_beta
cfg.tkmain =46;  % step size in hlselo_lqa
cfg.tol = 1e-5;
cfg.seed = 1;
cfg.verbose = true;
cfg.save_results = true;
cfg.out_dir = fullfile(pwd, 'results_hlselo');
end


function out = apply_defaults(in, defaults)
out = defaults;
if isempty(in)
    return;
end
f = fieldnames(in);
for k = 1:numel(f)
    out.(f{k}) = in.(f{k});
end
end


function [Z, X, T, C, Iota, Delta, R] = survival_data(n, Beta, mu, sigma, tau, iter)
% Z: n x p, X/T/C: n x 1, R: n x n
rng(iter, 'twister');
Z = mvnrnd(mu, sigma, n);

Beta2 = -Beta;
F = rand(n, 1);
lambda = exp(Z * Beta) + exp(Z * Beta2);
T = -log(1 - F) ./ lambda;
C = tau * rand(n, 1);
Delta = (T <= C);
X = min(T, C);

prop = exp(Z * Beta) ./ lambda;
Iota = double(rand(n, 1) < prop);

% Risk set R(j,i) = 1 if subject j at risk at time i
Xj = X(:);
Xi = X(:)';
cond1 = Xj >= Xi;
event_row = repmat((Delta(:) == 1 & Iota(:) ~= 1), 1, n);
cond2 = (Xj <= Xi) & event_row;
R = double(cond1 | cond2);
end


function W = Weight(X, T, C, Delta, n)
% Vectorized implementation of original weighting rule
G = C_KM(X, Delta);

Xi = X(:)';
Xj = X(:);
Ti = T(:)';
Ci = C(:)';

cond = Xi < Xj; % row j, col i
base = double(Ci >= bsxfun(@min, Ti, Xj));
ratio = G(:) ./ (G(:)' + 1e-6); % row j / col i

W = base;
W(cond) = base(cond) .* ratio(cond);

% keep original shape expectation explicit
if ~isequal(size(W), [n, n])
    error('Weight matrix size mismatch.');
end
end


function G = C_KM(X, Delta)
% Kaplan-Meier style estimate of censoring survival in original index order
n = numel(X);
cens = 1 - Delta;
[~, ord] = sort(X);
cens_sorted = cens(ord);

G_sorted = zeros(n, 1);
temp_G = 1;
for i = 1:n
    if cens_sorted(i) == 1 && (n - i ~= 0)
        temp_G = temp_G * (1 - 1 / (n - i + 1));
        G_sorted(i) = temp_G;
    elseif cens_sorted(i) == 0 && (n - i ~= 0)
        G_sorted(i) = temp_G;
    elseif cens_sorted(i) == 1 && (n - i == 0)
        G_sorted(i) = 0;
    else
        G_sorted(i) = temp_G;
    end
end

G = zeros(n, 1);
G(ord) = G_sorted;
end


function initial_beta = ini_beta(Z, r, Delta, Iota, R, W, tk)
[n, ~] = size(Z);
beta = zeros(size(Z, 2), 1);

k = 1;
err = false;
% tk = 4;
while k <= 1000 && ~err
    tem1 = (W .* R) .* exp(Z * beta);
    sum_tem1 = sum(tem1, 1);
    f2 = (Z' * tem1) ./ (sum_tem1 + 1e-12);
    f2 = f2';

    L_prime = -sum((Delta .* Iota) .* (Z - f2), 1)' / n;

    beta1 = beta - L_prime / tk;
    w = beta1 - beta;
    err = (norm(w, 2)^2 <= r * max(norm(beta, 2)^2, 1));
    beta = beta1;
    k = k + 1;
end

initial_beta = beta;
end


function [opt_beta, opt_theta, opt_delta] = hlselo_lqa(n, ini_beta, Z, r, Delta, Iota, R, W, theta, delta, tk)
[~, p] = size(Z);
opt_beta = zeros(p, 1);
beta1 = zeros(p, 1);
beta = ini_beta;
opt_BIC = inf;
n0 = sum(Delta);
lambda0 = log(max(n0, 2));

for h = 1:length(delta)
    for j = 1:length(theta)
        k = 1;
        err = false;
        % tk = 46;
        beta = ini_beta;

        while k <= 1000 && ~err
            tem1 = (W .* R) .* exp(Z * beta);
            sum_tem1 = sum(tem1, 1);
            f2 = (Z' * tem1) ./ (sum_tem1 + 1e-12);
            f2 = f2';

            L_prime = -sum((Delta .* Iota) .* (Z - f2), 1)' / n;
            beta_tilde = beta - L_prime / tk;

            idx_keep = find(abs(beta_tilde) > lambda0 * delta(h));
            idx_zero = setdiff(1:p, idx_keep);

            beta1(idx_zero) = 0;
            if ~isempty(idx_keep)
                pen_diag = lambda0 * (2 * abs(beta(idx_keep)) * theta(j) + delta(h)) .* ...
                    (1 ./ (1 + delta(h) * abs(beta(idx_keep)) + theta(j) * beta(idx_keep).^2).^2 + 1e-6);
                W1 = diag(pen_diag);
                u = eye(length(idx_keep)) + 2 * W1 / tk;
                beta1(idx_keep) = u \ beta_tilde(idx_keep);
            end

            w = beta1 - beta;
            err = (norm(w, 2)^2 <= r * max(norm(beta, 2)^2, 1));
            beta = beta1;
            k = k + 1;
        end

        beta2 = beta1;
        tem2 = sum((W .* R) .* exp(Z * beta2), 1);
        ell = -sum((Delta .* Iota) .* (Z * beta2 - log(tem2 + 1e-12))) / n;
        sel = (beta2 ~= 0);
        BIC = ell + sum(sel) * log(n) / n;

        if BIC <= opt_BIC
            opt_BIC = BIC;
            opt_theta = theta(j);
            opt_delta = delta(h);
            opt_beta = beta2;
        end
    end
end
end



function cov_diag = cov_PSH(index, Delta, Iota, Z, W, R, delta_h, theta, beta)
% Fine-Gray style sandwich covariance under LQA:
%   Var(beta_C) ~= A_C^{-1} B_C A_C^{-1}
% where A_C = -ell''_CC(beta_hat) + W_lambda,CC(beta_hat)
% and   B_C = sum_i U_i,C(beta_hat) U_i,C(beta_hat)'.
%
% Note: This uses IPCW-weighted score contributions U_i, but does not add
% an explicit KM(G)-influence augmentation term.
[n, ~] = size(Z);
k = numel(index);

eta = Z * beta;
exp_eta = exp(eta);
DI = Delta .* Iota;

% Sum-scale components
A_info = zeros(k, k);
U = zeros(k, n);

for i = 1:n
    wi = (W(:, i) .* R(:, i)) .* exp_eta;
    s0 = sum(wi) + 1e-12;
    s1 = Z' * wi;

    s2 = Z' * (Z .* wi);

    mu_i = s1 / s0;
    info_i = (s2 / s0) - (mu_i * mu_i');

    A_info = A_info + DI(i) * info_i(index, index);
    U(:, i) = DI(i) * (Z(i, index)' - mu_i(index));
end

n0 = sum(Delta);
lambda0 = log(max(n0, 2));
pen_diag = lambda0 * (2 * abs(beta(index)) * theta + delta_h) .* ...
    (1 ./ (1 + delta_h * abs(beta(index)) + theta * beta(index).^2).^2 + 1e-6);
W_lambda = diag(pen_diag);

A = A_info + W_lambda;
U_center = U - mean(U, 2);
B = U_center * U_center';

% Ridge-stabilize if nearly singular
ridge = 1e-8;
A_reg = A + ridge * eye(k);

V = (A_reg \ B) / A_reg';
cov_diag = max(real(diag(V)), 0);

if numel(cov_diag) ~= k
    cov_diag = zeros(k, 1);
end
end

