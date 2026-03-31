tic
format short
clc;clear;close;

%% === model parameters ===
n = 300; p = 10;
N = 200;  % simulation times
rho_set = [0.25 0.5];  % correlation parameter
mu = zeros(p,1); a = (1:p);
ama = bsxfun(@minus,a,a'); 
rho = rho_set(1);
sigma = rho.^(abs(ama));

%% === 'Censoring Rate' Setting ===
%-- Case A (Weak signal scenario) tau1-->censoring rate=25%
tau1 = 1.6;  tau2 = 0.75;    % rho = 0.25
% tau1 = 1.65;  tau2 = 0.78;   % rho = 0.5

%-- Case B (Strong signal scenario) tau1-->censoring rate=25%
% tau1 = 1.10;  tau2 = 0.45;   % rho = 0.25
% tau1 = 1.2;  tau2 = 0.5;     % rho = 0.5

%% === 'gamma' & 'delta' parameters ===
% %-- Case A (Weak signal scenario)
c_value = 0.15; h = 0.25;

% theta = 0.75*(1.13:0.0035:1.5)*1.95*10^4;  
theta = (0.60:-0.005:0.50)*1.25*10^4;
delta = (0.05:-0.0025:0.015);  %Corr=0.785  MSE=0.00803


%%-- Case B (Strong signal scenario)
% theta = (0.60:-0.005:0.50)*1.25*10^4; 
% theta = (1.05:-0.001:0.900)*1.095*10^3; 
% delta = 1.5;


%% === 'CD Algorithm's parameter ===
block_size = 20; 

%% === True Beta ===
Beta = 0.5*[-1.0;1.0;0;0;0;-1.0;zeros(p-6,1)];   %% Case A (Weak signal scenario)
% Beta = [-1.0;1.0;0;0;0;-1.0;zeros(p-6,1)];       %% Case B (Strong signal scenario)


%% === Defination of all terms ===
Censorrate = zeros(N,1);
index = find(Beta~=0);  % location to the true value

% Initialize arrays to store the results 
initial_beta = zeros(p,N);
hlselo = zeros(p,N); 

% Initialize arrays to store results for optimal tunning parameters
opt_theta = zeros(1,N); 
opt_delta = zeros(1,N);

% initial_beta setting
Q = 300;
initial_beta_Q = zeros(p,Q);

for j = 1: Q
    j
    rng(j*100)   % % Set random seed
    [Z,X,T,C,Iota,Delta,R] = survival_data(n,Beta,mu,sigma,tau1,j);  % CR = 45% -> tau2 && 25% -> tau1
    W = Weight(X,T,C,Delta,n);
    %% initial beta
    initial_beta_Q(:,j) = ini_beta(Z,1e-5,Delta,Iota,R,W);
end


%% === Conduct Monte Carlo simulations ===

for iter = 1:N
    iter 
    rng(iter)   % % Set random seed

   
   

    %% Generante the simulation data
    [Z,X,T,C,Iota,Delta,R] = survival_data(n,Beta,mu,sigma,tau1,iter);  % CR = 45% -> tau2 && 25% -> tau1
    Censorrate(iter) = 1-mean(Delta);
    W = Weight(X,T,C,Delta,n);

    %% initial beta
    % initial_beta(:,iter) = ini_beta(Z,1e-5,Delta,Iota,R,W);  
    % initial_beta(:,iter) = ini_beta_cd(Z,1e-5,Delta,Iota,R,W,theta,block_size);  % % initial beta - CD algorithm  

    % initial_beta_hselo (mean value of "initial_beta_Q ")
    mean_initial_beta = mean(initial_beta_Q,2);
    
    %% All variable selection methods
    % [hlselo(:,iter),opt_theta(iter),opt_delta(iter)] = hlselo_lqa(n,initial_beta(:,iter),Z,1e-5,Delta,Iota,R,W,theta/0.5,delta,c_value,h);   % % H-LSELO
    [hlselo(:,iter),opt_theta(iter),opt_delta(iter)] = hlselo_lqa(n, mean_initial_beta, Z,1e-5,Delta,Iota,R,W,theta/0.5,delta,c_value,h);   % % H-LSELO

    %[hlselo_cd(:,iter),opt_theta3(iter)] = lselo_CD(n,initial_beta(:,iter),Z,1e-5,Delta,Iota,R,W,theta2,block_size,p_size);

    %% MSE of each simulation
    se_hlselo(iter) = (hlselo(:,iter)-Beta)'*sigma*(hlselo(:,iter)-Beta);   % % mse_lselo
    
    %% covariance of each simulation
    cov_beta(:,iter) = cov_PSH(index,Delta,Iota,Z,W,R,opt_delta(iter),opt_theta(iter),hlselo(:,iter));  % covariance of coefficient

    
end


%% === Assessment Criteria ===
% -------------------------------------------------------------------------
% % BIAS;
% % SE;
% % MSE：Mean squre error;
% % ASE;
% % ASE_std;
% % CP_beta.

% -------------------------------------------------------------------------

[mean(initial_beta_Q(index,:),2) mean(hlselo(index,:),2)]

% -------------------------------------------------------------------------
BIAS = mean(hlselo(index,:),2)-Beta(index);

% -------------------------------------------------------------------------
SE = [std(hlselo(1,:)) std(hlselo(2,:)) std(hlselo(6,:))]';

% -------------------------------------------------------------------------
ASE = mean(sqrt(cov_beta),2);
ASE_std = [std(sqrt(cov_beta(1,:))) std(sqrt(cov_beta(2,:))) std(sqrt(cov_beta(3,:)))]';

% -------------------------------------------------------------------------
lower = hlselo(index,:) - 1.96 * sqrt(cov_beta);       
upper = hlselo(index,:) + 1.96 * sqrt(cov_beta);       

in_interval = (lower <= Beta(index)) & (upper >= Beta(index));  
CP_beta = mean(in_interval, 2);          

% CP_beta = mean( and( hlselo(index,:)-1.96*sqrt(cov_beta) <= repmat(Beta(index),1,N),...
%     hlselo(index,:)+1.96*sqrt(cov_beta) >= repmat(Beta(index),1,N) ),2);   % % Confidence intervals


%%-------------------------------------------------------------------------
corr = sum((all(hlselo(index,:))).*(1-any(hlselo(setdiff(1:1:p, index),:))))/N;
MSE = mean(se_hlselo);
N_plus = sum(sum(hlselo(setdiff(1:1:p, index),:)~=0))/N;
N_minus = sum(sum(hlselo(index,:)==0))/N;
Size = sum(sum(hlselo(:,:)~=0))/N;

%% === Print the simulation results ===

Criteria = [BIAS SE ASE ASE_std CP_beta];
method_name = ["\beta_{1}=-0.5"; "\beta_{2}=0.5"; "\beta_{6}=-0.5"];
Crit = [method_name Criteria];

methods = Crit(:, 1);  
numData = round(str2double(Crit(:, 2:end)), 4);  % Other columns: double, 4 decimal places
T = table(methods, numData(:,1), numData(:,2), numData(:,3), ...
          numData(:,4), numData(:,5), ...
          'VariableNames', {'Parameters','BIAS','SE','ASE','ASE_std','CP_beta'});

disp(T)


%
Criteria1 = [corr MSE N_plus N_minus Size];
method_name1 = "HSELO";
Crit1 = [method_name1 Criteria1];

methods1 = Crit1(:, 1);  
numData1 = round(str2double(Crit1(:, 2:end)), 4);  % Other columns: double, 4 decimal places
T1 = table(methods1, numData1(:,1), numData1(:,2), numData1(:,3), ...
          numData1(:,4), numData1(:,5), ...
          'VariableNames', {'Method','Pcorr','MSE','FPR','PNR','Size'});

disp(T1)

CR = mean(Censorrate)  % % Censoring rate

time = toc   % % Running time



%% ============================================================
%                  SUBFUNCTIONS
% ===============%==%==%==%==%==%==%===========================
%                 survival_data()
% ===============%==%==%==%==%==%==%===========================
function [Z,X,T,C,Iota,Delta,R] = survival_data(n,Beta,mu,sigma,tau,iter)
%Z:n*p
%X:n*1
%R:n*n
%p = length(Beta);
Z = mvnrnd(mu,sigma,n);
Beta2 = -Beta;

F = unifrnd(0,1,[n,1]);
lambda = exp(Z*Beta)+exp(Z*Beta2);
T = -log(1-F)./lambda;  %T:事件1和事件2总体下的生存时间
C = unifrnd(0,tau,[n,1]);
Delta = (T<=C);
X = min(T,C);

prop = exp(Z*Beta)./lambda; 
alphabet = [1 0];
Iota = zeros(n,1);
for i=1:n
    Iota(i,1) = randsrc(1,1,[alphabet; prop(i,1) 1-prop(i,1)]);
end

%**************************************************
% 历险集：R()
%**************************************************
R = zeros(n,n);
for i=1:n                        
    for j=1:n
        if X(j,1)>=X(i,1)||(X(j,1)<=X(i,1)&&Delta(j,1)==1&&Iota(j,1)~=1)
            R(j,i)=1;
        end
    end
end

end


%% ===================================================
%                 Weight()
% ============================================================
function W = Weight(X,T,C,Delta,n)
W = zeros(n,n);
% C_Delta = 1-Delta;
[G,id] = C_KM(X,Delta);

for i=1:n
    for j=1:n
        if X(i) < X(j)
            W(j,i) = (C(i)>=min(T(i),X(j)))*G(id(j))/(G(id(i))+1e-6);
        else
            W(j,i) = C(i)>=min(T(i),X(j));
        end

    end
end

end

%% ===================================================
%                 C_KM()
% ============================================================
function [G,id]=C_KM(X,Delta)
% [f,x] = ecdf(C,'Censoring',Delta); % Kaplan-Meier estimation for censored time 
% G = 1-f;
% ecdf(C,'censoring',Delta,'function','survivor');  % Survival function with censored time

[n,~] = size(X);
cens = 1-Delta;
[~,id] = sort(X);
cens = cens(id);   % Censoring indicator variable

temp_G=1; G = zeros(n,1);

for i=1:n 
    if cens(i,1)==1 && (n-i~=0)
        temp_G = temp_G*(1-1/(n-i+1));
        G(i,1) = temp_G;
    elseif cens(i,1)==0 && (n-i~=0)
        G(i,1) = temp_G;
    
    elseif cens(i,1)==1 && (n-i==0)
        G(i,1) = 0;
    else
        G(i,1) = temp_G;
    end
end

end


%% ===================================================
%                 ini_beta()
% ============================================================
function initial_beta = ini_beta(Z,r,Delta,Iota,R,W)
[n,p] = size(Z);
beta = zeros(p,1);
f2 = zeros(n,p);

k = 1; err = 0; tk = 4;
while k<=1000&&err==0
    % k
    for i=1:n
        tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
        %tem1 = R(:,i).*exp(Z*beta);
        f2(i,:) = sum(tem1.*Z)/sum(tem1);
    end    
    L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;
  
    beta1 =  beta - L_prime/tk;
    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

% beta2 = beta1; 
% f2 = zeros(n,p);
% 
% k = 1;err = 0; tk = 2;
% while k<=1000&&err==0
%     %k
%     for i=1:n
%         tem1 = (W(:,i).*R(:,i)).*exp(Z*beta2);
%         %tem1 = R(:,i).*exp(Z*beta);
%         f2(i,:) = sum(tem1.*Z)/sum(tem1);
%     end    
%     L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;
% 
%     beta3 =  beta2 - L_prime/tk;
%     w = beta3-beta2;
%     err = norm(w,2)^2 <= r*norm(beta2,2)^2;
%     beta2 = beta3;
%     k = k+1;
% end
% 
% initial_beta = beta3;

initial_beta = beta;

end


%% ===================================================
%                 ini_beta_cd()
% ============================================================
function initial_beta = ini_beta_cd(Z,r,Delta,Iota,R,W,theta,block_size)
[n,p] = size(Z);
beta = zeros(p,1);
f2 = zeros(n,p);

n0 = sum(Delta);
lambda0 = log(n0);
opt_BIC = 1e+10;


%% Block Coordinate Descent Parameter Settings
% block_size = 10;
num_blocks = ceil(p / block_size);
blocks = cell(num_blocks, 1);

%% Constructing Block Index
for b = 1:num_blocks
    start_idx = (b-1)*block_size + 1;
    end_idx = min(b*block_size, p);
    blocks{b} = start_idx:end_idx;
end


k = 1;err = 0; tk = 18;
while k<=1000&&err==0
    % k
    % for i=1:n
    %     tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
    %     %tem1 = R(:,i).*exp(Z*beta);
    %     f2(i,:) = sum(tem1.*Z)/sum(tem1);
    % end

    % revised on 17/03/2025
    tem1 = (W.*R).*exp(Z*beta);
    sum_tem1 = sum( tem1, 1 );
    f2 = (Z'*tem1)./sum_tem1;
    f2 = f2';

    L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

    %% Block coordinate descent
    if p<=10
        beta1 =  beta - L_prime/tk;
    else
        for iter = 1:50
            % 随机块遍历
            block_order = randperm(num_blocks);
            %%
            for b = block_order
                current_block = blocks{b};
                beta_block = beta(current_block);
                beta_new1 = beta_block - L_prime(current_block)/tk;

                beta1(current_block,1) = beta_new1;

            end

            delta = norm(beta - beta1,2)^2 / (norm(beta,2)^2 + 1e-8);
            if delta < 1e-6
                break;
            end
        end
    end

    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;

    % initial_beta = beta;

end

for j = 1:length(theta)

k=1;err=0; tk = 16;
while k<=1000 && err==0
    %k
    W1 = lambda0*theta(j)*diag(1./(1+theta(j)*beta.^2).^2);
    u = eye(p) + 2*W1/tk;

    % for i=1:n
    %     tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
    %     %tem1 = R(:,i).*exp(Z*beta);
    %     f2(i,:) = sum(tem1.*Z)/sum(tem1);
    % end

    % revised on 17/03/2025
    tem1 = (W.*R).*exp(Z*beta);
    sum_tem1 = sum( tem1, 1 );
    f2 = (Z'*tem1)./sum_tem1;
    f2 = f2';

    L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

    beta_tilde = beta - L_prime/tk;
    beta1 = u\beta_tilde;
    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

beta2 = beta.*(abs(beta)>=2*1e-4);

% for i=1:n
%     tem2(i,1) = sum((W(:,i).*R(:,i)).*exp(Z*beta2));
% end

% revised on 17/03/2025
tem2 = sum( (W.*R).*exp(Z*beta),1 );


ell = -sum((Delta.*Iota).*(Z*beta2-log(tem2)))/n;
sel = beta2~=0;
BIC = ell+sum(sel)*log(n)/n;

if BIC<=opt_BIC
    opt_BIC = BIC;
    opt_beta = beta2;
end

end

initial_beta = opt_beta;


end




%% ===================================================
%                 hlselo_lqa()
% ============================================================
function [opt_beta,opt_theta,opt_delta] = hlselo_lqa(n,ini_beta,Z,r,Delta,Iota,R,W,theta,delta,c,h)
[~,p] = size(Z);
% f2 = zeros(n,p);
% tem2 = zeros(n,1);
opt_beta = zeros(p,1);
beta1 = zeros(p,1);
beta = ini_beta;  % % initial value
% beta = zeros(p,1);
opt_BIC = 1e+10;
n0 = sum(Delta);
lambda0 = log(n0);

for h = 1:length(delta)

for j = 1:length(theta)

    k=1;err=0; tk = 46; % beta = ini_beta;
    while k<=1000 && err==0
        %k
        % tt = abs(beta)-c;
        % v_c = (delta*c+theta(j)*c^2)/( (1+delta*c+theta(j)*c^2)^2 + 1e-6);
        % d_c = (delta+2*c*theta(j))/( theta(j)*( (1/theta(j)+delta*c+theta(j)*c^2)^2 ) + 1e-6);
        % alpha = d_c*h/(1-v_c + 1e-6);
        % 
        % W1 = lambda0*( theta(j)+delta./(2*abs(beta) + 1e-6) ).*diag( (1./(1+delta*abs(beta)+theta(j)*beta.^2).^2 + 1e-6).*(abs(beta)<=c) ) +...
        %      lambda0*((c<abs(beta))&(abs(beta)<c+h)).*diag( (1-v_c)/(2*h*abs(beta) + 1e-6) ).*( -6*tt.^2+6*tt + alpha*(3*tt.^2-4*tt+1) );
        
        % W1 = lambda0*( theta(j)+delta./(2*abs(beta) + 1e-6) ).*diag( (1./(1+delta*abs(beta)+theta(j)*beta.^2).^2 + 1e-6));  


        % for i=1:n
        %     tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
        %     %tem1 = R(:,i).*exp(Z*beta);
        %     f2(i,:) = sum(tem1.*Z)/sum(tem1);
        % end

        % revised on 17/03/2025
        tem1 = (W.*R).*exp(Z*beta);
        sum_tem1 = sum( tem1, 1 );
        f2 = (Z'*tem1)./sum_tem1;
        f2 = f2';

        L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;


        beta_tilde = beta - L_prime/tk;
        index = find(abs(beta_tilde)>lambda0*delta(h));
        index1 = setdiff(1:p,index);
        beta1(index1) = zeros(length(index1),1);


        W1 = lambda0*( 2*abs(beta(index))*theta(j)+delta(h) ).*diag( (1./(1+delta(h)*abs(beta(index))+theta(j)*beta(index).^2).^2 + 1e-6));  
        % u = eye(p) + 2*W1/tk;
        u = eye(length(index)) + 2*W1/tk;

        beta1(index) = u\beta_tilde(index);

        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end
   
    % beta2 = beta.*(abs(beta)>=2*1e-4);
     beta2 = beta1;

    % for i=1:n
    %     tem2(i,1) = sum((W(:,i).*R(:,i)).*exp(Z*beta2));
    % end

    % revised on 17/03/2025
    tem2 = sum( (W.*R).*exp(Z*beta2),1 );

    ell = -sum((Delta.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;

    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_theta = theta(j);
        opt_delta = delta(h);
        opt_beta = beta2;
    end

end

end

end



%% ===================================================
%                 lselo_CD()
% ============================================================
function [opt_beta,opt_theta] = lselo_CD(n,ini_beta,Z,r,Delta,Iota,R,W,theta,block_size,p_size)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
opt_beta = zeros(p,1);
% beta = ini_beta;  % % Initial iteration value
% beta1 = zeros(p,1);
opt_BIC = 1e+10;
n0 = sum(Delta);
lambda0 = log(n0);


%% Block Coordinate Descent Parameter Settings
% block_size = 10;
num_blocks = ceil(p / block_size);
blocks = cell(num_blocks, 1);

%% Construct Block Index
for b = 1:num_blocks
    start_idx = (b-1)*block_size + 1;
    end_idx = min(b*block_size, p);
    blocks{b} = start_idx:end_idx;
end


for j = 1:length(theta)
    k=1; err=0; tk = 8; beta = ini_beta;
    while k<=1000 && err==0
        %k
        W1 = lambda0*theta(j)*diag(1./(1+theta(j)*beta.^2).^2);
        u = eye(p) + 2*W1/tk;

        % for i=1:n
        %     tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
        %     %tem1 = R(:,i).*exp(Z*beta);
        %     f2(i,:) = sum(tem1.*Z)/sum(tem1);
        % end

        % revised on 17/03/2025
        tem1 = (W.*R).*exp(Z*beta);
        sum_tem1 = sum( tem1, 1 );
        f2 = (Z'*tem1)./sum_tem1;
        f2 = f2';

        L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

        %% Block coordinate descent
        if p<=p_size
            beta_tilde = beta - L_prime/tk;
            % beta_tilde = beta - L_prime/(n*tk);
            beta1 = u\beta_tilde;
        else
            for iter = 1:50
                % Random block traversal
                block_order = randperm(num_blocks);
                %%
                for b = block_order
                    current_block = blocks{b};
                    beta_block = beta(current_block);

                    % Calculate W1
                    W12 = lambda0*theta(j)*diag(1./(1+theta(j)*beta_block.^2).^2);             
                    u1 = eye(length(current_block)) + 2*W12/tk;

                    beta_tilde = beta_block - L_prime(current_block)/tk;
                    beta_new1 = u1\beta_tilde;

                    beta1(current_block,1) = beta_new1;

                end

                delta = norm(beta - beta1,2)^2 / (norm(beta,2)^2 + 1e-8);
                if delta < 1e-6
                    break;
                end
            end
        end

        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;

        k = k+1;

    end
   
    beta2 = beta.*(abs(beta)>=2*1e-4);

    % for i=1:n
    %     tem2(i,1) = sum((W(:,i).*R(:,i)).*exp(Z*beta2));
    % end

    % revised on 17/03/2025
    tem2 = sum( (W.*R).*exp(Z*beta),1 );

    ell = -sum((Delta.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;

    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_theta = theta(j);
        opt_beta = beta2;
    end

end


end






%% ===================================================
%                 cov_PSH()
% ============================================================
%function cov_beta = cov_PSH(Z,r,Delta,Iota,R,W);

function cov_beta = cov_PSH(index,Delta,Iota, Z, W, R, delta_h, theta, beta)
[n,p] = size(Z);
tt = zeros(n,1);  % % n*1
Cel0 = cell(1,n); Cel1 = cell(1,n);
Cel2 = cell(1,n); Cel3 = cell(1,n);
Cel4 = cell(1,n);

eta = Z*beta;

% 计算分母
for i=1:n
    t = W(i,:)'.*exp(eta(i));
    tt(i) = sum(R(i,:)'.*t);
    Cel0{1,i} = Z(i,:)'*Z(i,:);
end

for i=1:n
    Cel0{1,i} = Z(i,:)'*Z(i,:);
end


% 计算分子
for j =1:n

    PP = (W(j,:).*R(j,:)).*exp(eta(j));

    for i=1:n
        KK = Cel0{1, i};
        Cel1{1,i} = PP(i)*KK;
        Cel2(1,i) = { ( ( (W(j,:).*R(j,:))'.*exp(eta(j)) ).*Z(j,:) )'*...
            ( ( (W(j,:).*R(j,:))'.*exp(eta(j)) ).*Z(j,:) ) };
    end

    M1 = zeros(p,p);
    M2 = zeros(p,p);

    for u=1:n
        M1 = M1 + Cel1{1,u};
        M2 = M2 + Cel2{1,u};
    end

    Cel3{1,j} = M1;
    Cel4{1,j} = M2;

end

DI = Delta.*Iota;
f1 = zeros(p,p);
f2 = zeros(p,p);

for i=1:n
  f1 = f1 + ( DI(i)/tt(i) )*Cel3{1,i};   
  f2 = f2 + ( DI(i)/tt(i)^2 )*Cel4{1,i}; 
end


% f1 = cellfun(@sum, Cel1, 'UniformOutput', false);   
% f2 = cellfun(@sum, Cel2, 'UniformOutput', false);  


L_primeprime = 2*(-f1+f2)/n;

% cov_beta = diag( inv( L_primeprime(index,index) ) );  

n0 = sum(Delta);
lambda0 = log(n0);

W1 = lambda0*( 2*abs(beta(index))*theta+delta_h ).*diag( (1./(1+delta_h*abs(beta(index))+theta*beta(index).^2).^2 + 1e-6));

term = L_primeprime(index,index)-W1;

cov_beta = diag( ( term\inv( L_primeprime(index,index)) )/term );  


end

Delete Sced_estimation.m (not needed)
