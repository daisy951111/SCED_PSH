tic
format short
clc; clear; close all;

%% === model parameters ===
n = 300;
p = 500;
N = 100;   % Monte Carlo times
rho_set = [0.25 0.5];
rho = rho_set(1);

mu = zeros(p,1);
a = 1:p;
ama = bsxfun(@minus, a, a');
sigma = rho.^(abs(ama));

%% === censoring setting ===
% Weak signal scenario
% rho=0.25 -> tau1=1.6 ; rho=0.5 -> tau1=1.65
tau1 = 1.6;

% Strong signal scenario examples
% tau1 = 1.10; % rho = 0.25
% tau1 = 1.20; % rho = 0.5

%% === theta & delta ===
theta = 0.95*(1.13:0.005:1.5)*1.95e4;
delta = 0.071;     % scalar or vector

%% === blockwise CD parameters ===
block_size = 10;
cd_inner_maxit = 50;

%% === Step size ===
tk_ini = 20;    % step size in ini_beta_cd
tk_main = 240;  % step size hlselo_CD

%% === true beta ===
Beta = 0.5*[-1.0;1.0;0;0;0;-1.0;zeros(p-6,1)];   % weak
% Beta = [-1.0;1.0;0;0;0;-1.0;zeros(p-6,1)];      % strong

%% === output control ===
do_cov = true;       % set true to compute covariance/CI
save_output = false;     % set true to save .mat/.csv
out_dir = fullfile(pwd, 'results_hlselo_cd_v2');

%% === containers ===
Censorrate = zeros(N,1);
index = find(Beta~=0);

initial_beta = zeros(p,N);
hlselo_cd = zeros(p,N);

se_hlselo = zeros(1,N);
cov_beta = zeros(numel(index),N);

opt_theta = zeros(1,N);
opt_delta = zeros(1,N);

% initial_beta setting
Q = 2;
initial_beta_R = zeros(p,Q);

%% === Monte Carlo ===
for iter = 1:N
    fprintf('iter = %d\n', iter);
    rng(iter, 'twister');

    for j = 1:Q
        % j
        rng(j*N*iter, 'twister')   % % Set random seed
        [Z_1,X,T,C,Iota,status,R] = survival_data(n,Beta,mu,sigma,tau1,j*N*iter);  % CR = 45% -> tau2 && 25% -> tau1
        W = Weight(X,T,C,status,n);

        %% initial beta
        initial_beta_R(:,j) = ini_beta_cd(Z_1, 1e-5, status, Iota, R, W, block_size, cd_inner_maxit, tk_ini);
    end

     % initial_beta_hselo (mean value of "initial_beta_Q ")
     initial_beta(:,iter) = mean(initial_beta_R,2);

    [Z,X,T,C,Iota,Delta,R] = survival_data(n,Beta,mu,sigma,tau1,iter);
    Censorrate(iter) = 1 - mean(Delta);
    W = Weight(X,T,C,Delta,n);

    % % initial beta (blockwise CD)
    % initial_beta(:,iter) = ini_beta_cd(Z, 1e-5, Delta, Iota, R, W, block_size, cd_inner_maxit);

    % HSELO + blockwise CD + warm starts along tuning grid
    [hlselo_cd(:,iter), opt_theta(iter), opt_delta(iter)] = hlselo_CD( ...
        n, initial_beta(:,iter), Z, 1e-5, Delta, Iota, R, W, theta, delta, block_size, cd_inner_maxit, tk_main);

    % MSE
    se_hlselo(iter) = (hlselo_cd(:,iter)-Beta)' * sigma * (hlselo_cd(:,iter)-Beta);

    % optional covariance
    if do_cov
        cov_beta(:,iter) = cov_PSH(index, Delta, Iota, Z, W, R, opt_delta(iter), opt_theta(iter), hlselo_cd(:,iter));
    end
end

%% === summary ===
disp([mean(initial_beta(index,:),2), mean(hlselo_cd(index,:),2)]);

Pcorr = mean( all(hlselo_cd(index,:)~=0,1) & all(hlselo_cd(setdiff(1:p,index),:)==0,1) );
MSE = mean(se_hlselo);
FPR = sum(sum(hlselo_cd(setdiff(1:p,index),:)~=0))/N;
PNR = sum(sum(hlselo_cd(index,:)==0))/N;
Size = sum(sum(hlselo_cd~=0))/N;

T1 = table("HSELO", Pcorr, MSE, FPR, PNR, Size, ...
    'VariableNames', {'Method','Pcorr','MSE','FPR','PNR','Size'});
disp(T1)

CR = mean(Censorrate);
fprintf('CR = %.4f\n', CR);
fprintf('time = %.2f sec\n', toc);

if do_cov
    BIAS = mean(hlselo_cd(index,:),2)-Beta(index);
    SE = std(hlselo_cd(index,:),0,2);
    ASE = mean(sqrt(max(cov_beta,0)),2);
    ASE_std = std(sqrt(max(cov_beta,0)),0,2);

    lower = hlselo_cd(index,:) - 1.96 * sqrt(max(cov_beta,0));
    upper = hlselo_cd(index,:) + 1.96 * sqrt(max(cov_beta,0));
    CP_beta = mean((lower <= Beta(index)) & (upper >= Beta(index)),2);

    Tcoef = table(string("beta_"+index), BIAS, SE, ASE, ASE_std, CP_beta, ...
        'VariableNames', {'Parameters','BIAS','SE','ASE','ASE_std','CP_beta'});
    disp(Tcoef)
end

if save_output
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    tag = sprintf('hlselo_cd_v2_n%d_p%d_N%d_rho%.2f_seed1', n, p, N, rho);
    save(fullfile(out_dir, [tag '.mat']), 'n','p','N','rho','tau1','theta','delta','block_size', ...
        'Beta','Censorrate','initial_beta','hlselo_cd','se_hlselo','opt_theta','opt_delta','cov_beta', ...
        'Pcorr','MSE','FPR','PNR','Size','CR','T1');
    writetable(T1, fullfile(out_dir, [tag '_model.csv']));
    if do_cov
        writetable(Tcoef, fullfile(out_dir, [tag '_coef.csv']));
    end
end

%% ============================================================
%                         SUBFUNCTIONS
%% ============================================================
function [Z,X,T,C,Iota,Delta,R] = survival_data(n,Beta,mu,sigma,tau,iter)
rng(iter, 'twister');
Z = mvnrnd(mu, sigma, n);
Beta2 = -Beta;

F = rand(n,1);
lambda = exp(Z*Beta) + exp(Z*Beta2);
T = -log(1-F)./lambda;
C = tau * rand(n,1);
Delta = (T<=C);
X = min(T,C);

prop = exp(Z*Beta)./lambda;
Iota = double(rand(n,1) < prop);

% Risk set R(j,i) = 1 if subject j at risk at X(i)
Xj = X(:);
Xi = X(:)';
cond1 = Xj >= Xi;
event_row = repmat((Delta(:)==1 & Iota(:)~=1), 1, n);
cond2 = (Xj <= Xi) & event_row;
R = double(cond1 | cond2);
end

function W = Weight(X,T,C,Delta,n)
G = C_KM(X,Delta);

Xi = X(:)';
Xj = X(:);
Ti = T(:)';
Ci = C(:)';

cond = Xi < Xj;
base = double(Ci >= bsxfun(@min, Ti, Xj));
ratio = G(:) ./ (G(:)' + 1e-6);

W = base;
W(cond) = base(cond) .* ratio(cond);

if ~isequal(size(W), [n,n])
    error('Weight matrix size mismatch.');
end
end

function G = C_KM(X,Delta)
n = numel(X);
cens = 1 - Delta;
[~,ord] = sort(X);
cens_sorted = cens(ord);

G_sorted = zeros(n,1);
temp_G = 1;
for i = 1:n
    if cens_sorted(i)==1 && (n-i~=0)
        temp_G = temp_G*(1-1/(n-i+1));
        G_sorted(i) = temp_G;
    elseif cens_sorted(i)==0 && (n-i~=0)
        G_sorted(i) = temp_G;
    elseif cens_sorted(i)==1 && (n-i==0)
        G_sorted(i) = 0;
    else
        G_sorted(i) = temp_G;
    end
end

G = zeros(n,1);
G(ord) = G_sorted;
end

function initial_beta = ini_beta_cd(Z,r,Delta,Iota,R,W,block_size,inner_maxit,tk)
[n,p] = size(Z);
beta = zeros(p,1);
beta1 = beta;

% tk = 20;
num_blocks = ceil(p/block_size);
blocks = make_blocks(p, block_size);

k = 1;
converged = false;
while k <= 1000 && ~converged
    % gradient at current beta
    L_prime = fg_score(beta, Z, Delta, Iota, R, W, n);

    beta1 = beta;
    if p <= block_size
        beta1 = beta - L_prime/tk;
    else
        for it = 1:inner_maxit
            beta_old_inner = beta1;
            % recompute score after each full sweep for better CD stability
            L_prime_inner = fg_score(beta1, Z, Delta, Iota, R, W, n);
            for b = randperm(num_blocks)
                idx = blocks{b};
                beta1(idx) = beta1(idx) - L_prime_inner(idx)/tk;
            end
            if norm(beta1-beta_old_inner,2)^2 <= 1e-6*max(norm(beta_old_inner,2)^2,1)
                break;
            end
        end
    end

    w = beta1 - beta;
    converged = norm(w,2)^2 <= r*max(norm(beta,2)^2,1);
    beta = beta1;
    k = k + 1;
end

initial_beta = beta1;
end

function [opt_beta,opt_theta,opt_delta] = hlselo_CD(n,ini_beta,Z,r,Delta,Iota,R,W,theta,delta,block_size,inner_maxit,tk)
[~,p] = size(Z);

opt_beta = zeros(p,1);
opt_BIC = inf;

n0 = sum(Delta);
lambda0 = log(max(n0,2));

num_blocks = ceil(p/block_size);
blocks = make_blocks(p, block_size);

% warm start between tuning points
beta_warm = ini_beta;

for h = 1:numel(delta)
    for j = 1:numel(theta)
        beta = beta_warm;
        beta1 = beta;
        % tk = 240;
        k = 1;
        converged = false;

        while k <= 1000 && ~converged
            % outer score at current beta
            L_prime = fg_score(beta, Z, Delta, Iota, R, W, n);

            if p <= block_size
                beta_tilde = beta - L_prime / tk;
                idx_keep = find(abs(beta_tilde) > lambda0*delta(h));
                idx_zero = setdiff(1:p, idx_keep);
                beta1 = beta;
                beta1(idx_zero) = 0;

                if ~isempty(idx_keep)
                    pen = lambda0 * (2 * abs(beta(idx_keep)) * theta(j) + delta(h)) .* ...
                        (1 ./ (1 + delta(h)*abs(beta(idx_keep)) + theta(j)*beta(idx_keep).^2).^2 + 1e-6);
                    U = eye(numel(idx_keep)) + 2*diag(pen)/tk;
                    beta1(idx_keep) = U \ beta_tilde(idx_keep);
                end
            else
                beta1 = beta;
                for it = 1:inner_maxit
                    beta_old_inner = beta1;
                    % refresh score at latest beta1 for better block updates
                    L_prime_inner = fg_score(beta1, Z, Delta, Iota, R, W, n);

                    for b = randperm(num_blocks)
                        idx = blocks{b};
                        bcur = beta1(idx);
                        btilde = bcur - L_prime_inner(idx)/tk;

                        keep = find(abs(btilde) > lambda0*delta(h));
                        bnew = zeros(numel(idx),1);

                        if ~isempty(keep)
                            bkeep = bcur(keep);
                            pen = lambda0 * (2 * abs(bkeep) * theta(j) + delta(h)) .* ...
                                (1 ./ (1 + delta(h)*abs(bkeep) + theta(j)*bkeep.^2).^2 + 1e-6);
                            U = eye(numel(keep)) + 2*diag(pen)/tk;
                            bnew(keep) = U \ btilde(keep);
                        end

                        beta1(idx) = bnew;
                    end

                    if norm(beta1-beta_old_inner,2)^2 <= 1e-6*max(norm(beta_old_inner,2)^2,1)
                        break;
                    end
                end
            end

            w = beta1 - beta;
            converged = norm(w,2)^2 <= r*max(norm(beta,2)^2,1);
            beta = beta1;
            k = k + 1;
        end

        beta2 = beta1;
        tem2 = sum((W.*R).*exp(Z*beta2),1) + 1e-12;
        ell = -sum((Delta.*Iota).*(Z*beta2 - log(tem2)'))/n;
        BIC = ell + nnz(beta2)*log(n)/n;

        if BIC <= opt_BIC
            opt_BIC = BIC;
            opt_theta = theta(j);
            opt_delta = delta(h);
            opt_beta = beta2;
        end

        beta_warm = beta2;
    end
end
end

function g = fg_score(beta, Z, Delta, Iota, R, W, n)
tem1 = (W.*R).*exp(Z*beta);
sum_tem1 = sum(tem1,1) + 1e-12;
f2 = (Z' * tem1) ./ sum_tem1;
f2 = f2';
g = -sum((Delta.*Iota).*(Z-f2),1)'/n;
end

function blocks = make_blocks(p, block_size)
num_blocks = ceil(p/block_size);
blocks = cell(num_blocks,1);
for b = 1:num_blocks
    s = (b-1)*block_size + 1;
    e = min(b*block_size, p);
    blocks{b} = s:e;
end
end

function cov_beta = cov_PSH(index,Delta,Iota,Z,W,R,delta_h,theta,beta)
% LQA + Fine-Gray style sandwich on selected set
[n,~] = size(Z);
k = numel(index);

eta = Z*beta;
exp_eta = exp(eta);
DI = Delta.*Iota;

A_info = zeros(k,k);
U = zeros(k,n);

for i = 1:n
    wi = (W(:,i).*R(:,i)).*exp_eta;
    s0 = sum(wi) + 1e-12;
    s1 = Z' * wi;
    s2 = Z' * (Z .* wi);

    mu_i = s1/s0;
    info_i = (s2/s0) - (mu_i*mu_i');

    A_info = A_info + DI(i)*info_i(index,index);
    U(:,i) = DI(i) * (Z(i,index)' - mu_i(index));
end

lambda0 = log(max(sum(Delta),2));
pen_diag = lambda0*(2*abs(beta(index))*theta + delta_h) .* ...
    (1./(1 + delta_h*abs(beta(index)) + theta*beta(index).^2).^2 + 1e-6);

A = A_info + diag(pen_diag);
U_center = U - mean(U,2);
B = U_center * U_center';

ridge = 1e-8;
A_reg = A + ridge*eye(k);
V = (A_reg \ B) / A_reg';

cov_beta = max(real(diag(V)),0);
end