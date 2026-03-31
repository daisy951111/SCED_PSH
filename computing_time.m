tic
format short
clc;clear;close;

%% === 'model parameters' ===
n1 = 300; n2 = 500; p1 = 10; p2 = 40;
N = 100;  % simulation times
rho_set = [0.25 0.5];  % correlation parameter
mu1 = zeros(p1,1); a1 = (1:p1);
mu2 = zeros(p2,1); a2 = (1:p2);
ama1 = bsxfun(@minus,a1,a1'); 
ama2 = bsxfun(@minus,a2,a2'); 
rho = rho_set(1);
sigma1 = rho.^(abs(ama1));
sigma2 = rho.^(abs(ama2));

%% === 'Censoring Rate' Setting ===
%-- Case A (Weak signal scenario) tau1-->censoring rate=25%
tau1 = 1.6;  tau2 = 0.75;    % rho = 0.25
% tau1 = 1.65;  tau2 = 0.78;   % rho = 0.5

%-- Case B (Strong signal scenario) tau1-->censoring rate=25%
% tau1 = 1.10;  tau2 = 0.45;   % rho = 0.25
% tau1 = 1.2;  tau2 = 0.5;     % rho = 0.5

%% === 'gamma' & 'delta' parameters ===
% %-- Case A (Weak signal scenario)

% theta = 0.75*(1.13:0.0035:1.5)*1.95*10^4;  
theta = (0.60:-0.0095:0.50)*1.25*10^4;


delta = (0.055:-0.009:0.015);  %Corr=0.785  MSE=0.00803


%%-- Case B (Strong signal scenario)
% theta = (0.60:-0.005:0.50)*1.25*10^4; 
% theta = (1.05:-0.001:0.900)*1.095*10^3; 
% delta = 1.5;


%% === 'lambda' parameter ===
lambda_selo = 0.15*(0.001).^((1:100)/100);
lambda_MCP = 0.15*(0.001).^((1:100)/100);
lambda_scad = 0.15*(0.001).^((1:100)/100);
lambda_Alasso = 0.15*(0.001).^((1:100)/100);
lambda_lasso = 0.15*(0.001).^((1:100)/100);

%% === 'theta' parameter for gselo ===
% alpha = (1.35:-0.005:0.900)*10^2; 
alpha = (1.35:-0.002:0.900)*10^2; 

%% === 'a' parameter for MIC ===
% a = 90;
a = (70:1:200);

%% === 'CD Algorithm's parameter ===
block_size = 50; 

%% === True Beta ===
Beta1 = 0.5*[-1.0;1.0;0;0;0;-1.0;zeros(p1-6,1)];   %% Case A (Weak signal scenario)
Beta2 = 0.5*[-1.0;1.0;0;0;0;-1.0;zeros(p2-6,1)];
% Beta = [-1.0;1.0;0;0;0;-1.0;zeros(p-6,1)];       %% Case B (Strong signal scenario)


%% === Defination of all terms ===
Censorrate = zeros(N,1);
% index = find(Beta~=0);  % location to the true value

% Initialize arrays to store the results 
% initial_beta = zeros(p,N);
% hlselo_cd = zeros(p,N); 


% Initialize arrays to store results for optimal tunning parameters
opt_theta = zeros(1,N); 
opt_delta = zeros(1,N);


%% === Conduct Monte Carlo simulations ===

for iter = 1:N

    iter 
    rng(iter)   % % Set random seed

    %% Generante the simulation data
    [Z,X,T,C,Iota,Delta,R] = survival_data(n1,Beta1,mu1,sigma1,tau1,iter);  % CR = 45% -> tau2 && 25% -> tau1
    Censorrate(iter) = 1-mean(Delta);
    W = Weight(X,T,C,Delta,n1);

    %% initial beta
  
    initial_beta1(:,iter) = ini_beta_cd(Z,1e-5,Delta,Iota,R,W,block_size);  % % initial beta - CD algorithm  
    
    %% All variable selection methods
    %--SCED methods
    tic;
    [hlselo_cd(:,iter),opt_theta(iter),opt_delta(iter)] = hlselo_CD(n1,initial_beta1(:,iter),Z,1e-5,Delta,Iota,R,W,theta/0.5,delta,block_size);
    times1(iter, 1) = toc;

    % se_hlselo(iter) = (hlselo_cd(:,iter)-Beta)'*sigma*(hlselo_cd(:,iter)-Beta);   % % mse_lselo

    %--Other methods
    tic;
    [gselo(:,iter),opt_alpha(iter)] = gselo_lqa(n1,initial_beta1(:,iter),Z,1e-5,Delta,Iota,R,W,alpha);   % % GSELO
    times1(iter, 2) = toc;

    %--MIC
    tic;
    [MIC(:,iter),opt_a(iter)] = MIC_lqa(n1,initial_beta1(:,iter),Z,1e-5,Delta,Iota,R,W,a);    
    times1(iter, 3) = toc;

    %--SELO
    tic;
    [selo(:,iter),opt_lambda_selo(iter)] = SELO_lqa(n1,initial_beta1(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_selo);     
    times1(iter, 4) = toc;

    %--SCAD
    tic;
    [scad(:,iter),opt_lambda_scad(iter)] = scad_lqa(n1,initial_beta1(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_scad);   
    times1(iter, 5) = toc;

    %--MCP  
    tic;
    [MCP(:,iter),opt_lambda_MCP(iter)] = MCP_lqa(n1,initial_beta1(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_MCP);       
    times1(iter, 6) = toc;

    %--Alasso
    tic;
    [Alasso(:,iter),opt_lambda_Alasso(iter)] = Alasso_lqa(n1,initial_beta1(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_Alasso);  
    times1(iter, 7) = toc;

    %--Lasso
    tic;
    [lasso(:,iter),opt_lambda_lasso(iter)] = lasso_lqa(n1,initial_beta1(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_lasso);      
    times1(iter, 8) = toc;



    %% Generante the simulation data
    [Z,X,T,C,Iota,Delta,R] = survival_data(n1,Beta2,mu2,sigma2,tau1,iter);  % CR = 45% -> tau2 && 25% -> tau1
    Censorrate(iter) = 1-mean(Delta);
    W = Weight(X,T,C,Delta,n1);

    %% initial beta
  
    initial_beta2(:,iter) = ini_beta_cd(Z,1e-5,Delta,Iota,R,W,block_size);  % % initial beta - CD algorithm  
    
    %% All variable selection methods
    %--SCED methods
    tic;
    [hlselo_cd2(:,iter),opt_theta2(iter),opt_delta2(iter)] = hlselo_CD(n1,initial_beta2(:,iter),Z,1e-5,Delta,Iota,R,W,theta/0.5,delta,block_size);
    times2(iter, 1) = toc;

    % se_hlselo(iter) = (hlselo_cd(:,iter)-Beta)'*sigma*(hlselo_cd(:,iter)-Beta);   % % mse_lselo

    %--Other methods
    tic;
    [gselo2(:,iter),opt_alpha2(iter)] = gselo_lqa(n1,initial_beta2(:,iter),Z,1e-5,Delta,Iota,R,W,alpha);   % % GSELO
    times2(iter, 2) = toc;

    %--MIC
    tic;
    [MIC2(:,iter),opt_a2(iter)] = MIC_lqa(n1,initial_beta2(:,iter),Z,1e-5,Delta,Iota,R,W,a);    
    times2(iter, 3) = toc;

    %--SELO
    tic;
    [selo2(:,iter),opt_lambda_selo2(iter)] = SELO_lqa(n1,initial_beta2(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_selo);     
    times2(iter, 4) = toc;

    %--SCAD
    tic;
    [scad2(:,iter),opt_lambda_scad2(iter)] = scad_lqa(n1,initial_beta2(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_scad);   
    times2(iter, 5) = toc;

    %--MCP  
    tic;
    [MCP2(:,iter),opt_lambda_MCP(iter)] = MCP_lqa(n1,initial_beta2(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_MCP);       
    times2(iter, 6) = toc;

    %--Alasso
    tic;
    [Alasso2(:,iter),opt_lambda_Alasso2(iter)] = Alasso_lqa(n1,initial_beta2(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_Alasso);  
    times2(iter, 7) = toc;

    %--Lasso
    tic;
    [lasso2(:,iter),opt_lambda_lasso2(iter)] = lasso_lqa(n1,initial_beta2(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_lasso);      
    times2(iter, 8) = toc;



    %% Generante the simulation data
    [Z,X,T,C,Iota,Delta,R] = survival_data(n2,Beta1,mu1,sigma1,tau1,iter);  % CR = 45% -> tau2 && 25% -> tau1
    Censorrate(iter) = 1-mean(Delta);
    W = Weight(X,T,C,Delta,n2);

    %% initial beta
  
    initial_beta3(:,iter) = ini_beta_cd(Z,1e-5,Delta,Iota,R,W,block_size);  % % initial beta - CD algorithm  
    
    %% All variable selection methods
    %--SCED methods
    tic;
    [hlselo_cd3(:,iter),opt_theta3(iter),opt_delta3(iter)] = hlselo_CD(n2,initial_beta3(:,iter),Z,1e-5,Delta,Iota,R,W,theta/0.5,delta,block_size);
    times3(iter, 1) = toc;

    % se_hlselo(iter) = (hlselo_cd(:,iter)-Beta)'*sigma*(hlselo_cd(:,iter)-Beta);   % % mse_lselo

    %--Other methods
    tic;
    [gselo3(:,iter),opt_alpha3(iter)] = gselo_lqa(n2,initial_beta3(:,iter),Z,1e-5,Delta,Iota,R,W,alpha);   % % GSELO
    times3(iter, 2) = toc;

    %--MIC
    tic;
    [MIC3(:,iter),opt_a3(iter)] = MIC_lqa(n2,initial_beta3(:,iter),Z,1e-5,Delta,Iota,R,W,a);    
    times3(iter, 3) = toc;

    %--SELO
    tic;
    [selo3(:,iter),opt_lambda_selo3(iter)] = SELO_lqa(n2,initial_beta3(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_selo);     
    times3(iter, 4) = toc;

    %--SCAD
    tic;
    [scad3(:,iter),opt_lambda_scad3(iter)] = scad_lqa(n2,initial_beta3(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_scad);   
    times3(iter, 5) = toc;

    %--MCP  
    tic;
    [MCP3(:,iter),opt_lambda_MCP3(iter)] = MCP_lqa(n2,initial_beta3(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_MCP);       
    times3(iter, 6) = toc;

    %--Alasso
    tic;
    [Alasso3(:,iter),opt_lambda_Alasso3(iter)] = Alasso_lqa(n2,initial_beta3(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_Alasso);  
    times3(iter, 7) = toc;

    %--Lasso
    tic;
    [lasso3(:,iter),opt_lambda_lasso3(iter)] = lasso_lqa(n2,initial_beta3(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_lasso);      
    times3(iter, 8) = toc;



    %% Generante the simulation data
    [Z,X,T,C,Iota,Delta,R] = survival_data(n2,Beta2,mu2,sigma2,tau1,iter);  % CR = 45% -> tau2 && 25% -> tau1
    Censorrate(iter) = 1-mean(Delta);
    W = Weight(X,T,C,Delta,n2);

    %% initial beta
  
    initial_beta4(:,iter) = ini_beta_cd(Z,1e-5,Delta,Iota,R,W,block_size);  % % initial beta - CD algorithm  
    
    %% All variable selection methods
    %--SCED methods
    tic;
    [hlselo_cd4(:,iter),opt_theta4(iter),opt_delta4(iter)] = hlselo_CD(n2,initial_beta4(:,iter),Z,1e-5,Delta,Iota,R,W,theta/0.5,delta,block_size);
    times4(iter, 1) = toc;

    % se_hlselo(iter) = (hlselo_cd(:,iter)-Beta)'*sigma*(hlselo_cd(:,iter)-Beta);   % % mse_lselo

    %--Other methods
    tic;
    [gselo4(:,iter),opt_alpha4(iter)] = gselo_lqa(n2,initial_beta4(:,iter),Z,1e-5,Delta,Iota,R,W,alpha);   % % GSELO
    times4(iter, 2) = toc;

    %--MIC
    tic;
    [MIC4(:,iter),opt_a4(iter)] = MIC_lqa(n2,initial_beta4(:,iter),Z,1e-5,Delta,Iota,R,W,a);    
    times4(iter, 3) = toc;

    %--SELO
    tic;
    [selo4(:,iter),opt_lambda_selo4(iter)] = SELO_lqa(n2,initial_beta4(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_selo);     
    times4(iter, 4) = toc;

    %--SCAD
    tic;
    [scad4(:,iter),opt_lambda_scad4(iter)] = scad_lqa(n2,initial_beta4(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_scad);   
    times4(iter, 5) = toc;

    %--MCP  
    tic;
    [MCP4(:,iter),opt_lambda_MCP4(iter)] = MCP_lqa(n2,initial_beta4(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_MCP);       
    times4(iter, 6) = toc;

    %--Alasso
    tic;
    [Alasso4(:,iter),opt_lambda_Alasso4(iter)] = Alasso_lqa(n2,initial_beta4(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_Alasso);  
    times4(iter, 7) = toc;

    %--Lasso
    tic;
    [lasso4(:,iter),opt_lambda_lasso4(iter)] = lasso_lqa(n2,initial_beta4(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_lasso);      
    times4(iter, 8) = toc;

    
end


%% === Assessment Criteria ===
%%-------------------------------------------------------------------------
% corr = sum((all(hlselo_cd(index,:))).*(1-any(hlselo_cd(setdiff(1:1:p, index),:))))/N;
% MSE = mean(se_hlselo);
% N_plus = sum(sum(hlselo_cd(setdiff(1:1:p, index),:)~=0))/N;
% N_minus = sum(sum(hlselo_cd(index,:)==0))/N;
% Size = sum(sum(hlselo_cd(:,:)~=0))/N;
% 
% %
% Criteria1 = [corr MSE N_plus N_minus Size];
% method_name1 = "HSELO";
% Crit1 = [method_name1 Criteria1];
% 
% methods1 = Crit1(:, 1);  
% numData1 = round(str2double(Crit1(:, 2:end)), 4);  % Other columns: double, 4 decimal places
% T1 = table(methods1, numData1(:,1), numData1(:,2), numData1(:,3), ...
%           numData1(:,4), numData1(:,5), ...
%           'VariableNames', {'Method','Pcorr','MSE','FPR','PNR','Size'});
% 
% disp(T1)

CR = mean(Censorrate)  % % Censoring rate

time_all = toc   % % Running time


% --- Data Analysis and Visualization ---

num_methods = length(times1(1,:));            % Number of methods to compare

% 1. Calculate Statistics
mean_times = mean(times1);       % Mean time
std_times = std(times1);         % Standard deviation
median_times = median(times1);   % Median time (Recommended for comparison as it is robust to outliers)

% 2. Print Results
method_names = {'SCED', 'GSELO', 'MIC', 'SELO', 'SCAD', 'MCP', 'ALASSO', 'LASSO'};
fprintf('\n=== Computation Time Comparison (Unit: Seconds) ===\n');
fprintf('%-12s | %-12s | %-12s | %-12s\n', 'Method', 'Mean', 'Median', 'Std Dev');
fprintf('------------------------------------------------------------\n');

for i = 1:num_methods
    fprintf('%-12s | %-12.6f | %-12.6f | %-12.6f\n', ...
        method_names{i}, mean_times(i), median_times(i), std_times(i));
end

% % 3. Visualization (Boxplot best shows distribution and outliers)
% figure;
% boxplot(times1, 'Labels', method_names);
% title('Distribution of Computation Times for 7 Methods over 200 Monte Carlo Simulations');
% ylabel('Time (seconds)');
% grid on;
% 
% % Alternatively, plot a bar chart of mean times
% figure;
% bar(mean_times);
% set(gca, 'XTickLabel', method_names);
% title('Average Computation Time Comparison');
% ylabel('Time (seconds)');
% grid on;


% --- Combined Visualization: 2x2 Grid Layout ---

% 1. Data Preparation
% Define the dataset cell array and corresponding titles
all_datasets = {times1, times2, times3, times4}; 
plot_titles = {'(a) n=300; p=10', '(b) n=300; p=40', '(c) n=500; p=10', '(d) n=500; p=40'};
% method_names = {'Method 1', 'Method 2', 'Method 3', 'Method 4', 'Method 5', 'Method 6', 'Method 7', 'Method 8'};
method_names = {'SCED', 'GSELO', 'MIC', 'SELO', 'SCAD', 'MCP', 'ALASSO', 'LASSO'};

% 1. 瀹氫箟鐗╃悊灏哄 (鑻卞)
% 鍙屾爮瀹藉害锛?7.0 鑻卞
figWidth = 8.5;
% 楂樺害锛氭牴鎹? 2x2 姣斾緥浼扮畻锛岀◢鍚庣敱 exportgraphics 鑷姩寰皟鎴栧浐瀹?
figHeight = 9.0; 

% 2. 鍒涘缓鍥剧獥
fig = figure('Color', 'w', ...
    'Units', 'inches', ...
    'Position', [1, 1, figWidth, figHeight], ...
    'Visible', 'on');


% 2. Create Figure with Tiled Layout (Requires MATLAB R2019b or newer)
% fig  = figure('Color', 'w', 'Position', [100, 100, 1000, 800]);
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');


% Add a global title for the entire figure
% title(t, 'Distribution of Computation Times Across 4 Simulation Sets', ...
%     'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% title(t,'Distribution of Computation Times for 8 Methods over 100 Monte Carlo Simulations', ...
%     'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% 3. Loop to Plot Each Dataset
for i = 1:4
    nexttile; % Move to the next tile in the 2x2 grid
    
    % --- Core Plotting Command (Based on your original code) ---
    % Original: boxplot(times1, 'Labels', method_names);
    % Adapted for loop:

    % colors = [
    % 0.12, 0.47, 0.71; % 钃?
    % 0.85, 0.33, 0.10; % 姗?
    % 0.00, 0.60, 0.40; % 缁?
    % 0.90, 0.20, 0.20; % 绾?
    % 0.50, 0.00, 0.50; % 绱?
    % 0.60, 0.40, 0.20; % 妫?
    % 0.90, 0.90, 0.90  % 鐏? (鐢ㄤ簬澶囩敤)
    % ];
    % 
    % 
    % boxplot(all_datasets{i}, 'Labels', method_names, 'Colors', colors, 'BoxStyle', 'filled');


    boxplot(all_datasets{i}, 'Labels', method_names);
    
    % Add specific title for this subplot
    title(plot_titles{i}, 'FontSize', 11, 'FontWeight', 'bold');
    
    % Basic Formatting
    ax = gca;
    ax.FontName = 'Arial';
    ax.FontSize = 9;
    ax.LineWidth = 0.8;
    ax.YGrid = 'off';      % Enable grid for better readability
    ax.GridAlpha = 0.3;   % Make grid lines faint
    
    % Axis Labels Management
    % Only show Y-label on the left column (plots 1 and 3) to avoid repetition
    if i == 1 || i == 3
        ylabel('Time (seconds)', 'FontSize', 10, 'FontWeight', 'bold');
    else
        ylabel('');
    end
    
    % Rotate X-axis labels if they overlap (optional, usually good for 7 methods)
    xtickangle(45); 
end


fprintf('2x2 Grid plot generated successfully.\n');

print(fig, '-dpdf', '-r600', 'comparison_4sets.eps');

 % --- 鏇挎崲鍘熸潵鐨? print 鍛戒护 ---
filename = 'comparison_4sets.pdf';

% 'Padding', 'tight' 鏄牳蹇冿細瀹冧細绱ц创鐫?浣犵殑鍧愭爣杞存爣绛捐鍓紝涓嶇暀浠讳綍搴曢儴鐧借竟
exportgraphics(fig, filename, 'ContentType', 'vector', 'Padding', 'tight');



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
T = -log(1-F)./lambda;  %T:浜嬩欢1鍜屼簨浠?2鎬讳綋涓嬬殑鐢熷瓨鏃堕棿
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
% 鍘嗛櫓闆嗭細R()
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
%                 ini_beta_cd()
% ============================================================
function initial_beta = ini_beta_cd(Z,r,Delta,Iota,R,W,block_size)
[n,p] = size(Z);
beta = zeros(p,1);
% f2 = zeros(n,p);

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


k = 1;err = 0; tk = 10;
while k<=1000&&err==0
    % k
    
    tem1 = (W.*R).*exp(Z*beta);
    sum_tem1 = sum( tem1, 1 );
    f2 = (Z'*tem1)./sum_tem1;
    f2 = f2';

    L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

    %% Block coordinate descent
    if p <= block_size
        beta1 =  beta - L_prime/tk;
    else
        for iter = 1:50
            % Random block traversal
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

end

initial_beta = beta1;


end





%% ===================================================
%                 hlselo__CD()
% ============================================================
function [opt_beta,opt_theta,opt_delta] = hlselo_CD(n,ini_beta,Z,r,Delta,Iota,R,W,theta,delta,block_size)

[~,p] = size(Z);
opt_beta = zeros(p,1);
% beta = ini_beta;  % % Initial iteration value
beta1 = zeros(p,1);
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

% beta_new1 = zeros(block_size,1);

for h = 1:length(delta)
    for j = 1:length(theta)
        k=1; err=0; tk = 80; beta = ini_beta;
        while k<=1000 && err==0
            %k

            tem1 = (W.*R).*exp(Z*beta);
            sum_tem1 = sum( tem1, 1 );
            f2 = (Z'*tem1)./sum_tem1;
            f2 = f2';

            L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

            %% Block coordinate descent
            if p <= block_size

                beta_tilde = beta - L_prime / tk;
                index = find(abs(beta_tilde) > lambda0*delta(h));
                index1 = setdiff(1:p,index);
                beta1(index1) = zeros(length(index1),1);

                W1 = lambda0 * (2 * abs(beta(index)) * theta(j) + delta(h)) .* ...
                    diag(1 ./ (1 + delta(h)*abs(beta(index)) + theta(j)*beta(index).^2).^2 + 1e-6);
                u = eye(length(index)) + 2 * W1 / tk;
                beta1(index) = u \ beta_tilde(index);

            else
                for iter = 1:50
                    % Random block traversal
                    block_order = randperm(num_blocks);
                    %%
                    for b = block_order
                        current_block = blocks{b};
                        block_len = length(current_block);  % 鍔ㄦ?佽幏鍙栧綋鍓嶅潡闀垮害
                        beta_block = beta(current_block);

                        % 鍔ㄦ?佸垱寤? beta_new1锛屽ぇ灏忓尮閰嶅綋鍓嶅潡
                        beta_new1 = zeros(block_len, 1);

                        beta_tilde_block = beta_block - L_prime(current_block) / tk;
                        index_block = find(abs(beta_tilde_block) > lambda0*delta(h));
                        index1_block = setdiff(1:length(beta_block),index_block);
                        beta_new1(index1_block) = zeros(length(index1_block),1);

                        % Calculate W2
                        W2 = lambda0 * (2 * abs(beta_block(index_block)) * theta(j) + delta(h)) .* ...
                            diag(1 ./ (1 + delta(h)*abs(beta_block(index_block)) + theta(j)*beta_block(index_block).^2).^2 + 1e-6);
                        u1 = eye(length(index_block)) + 2*W2/tk;

                        beta_new1(index_block) = u1 \ beta_tilde_block(index_block);
                        beta1(current_block,1) = beta_new1;

                    end

                    gamma1 = norm(beta - beta1,2)^2 / (norm(beta,2)^2 + 1e-8);

                    if gamma1 < 1e-6
                        break;
                    end
                end
            end

            w = beta1-beta;
            err = norm(w,2)^2 <= r*norm(beta,2)^2;
            beta = beta1;

            k = k+1;

        end

        % beta2 = beta.*(abs(beta)>=2*1e-4);
        beta2 = beta1;

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
%                 gselo_lqa()
% ============================================================
function [opt_beta,opt_theta] = gselo_lqa(n,ini_beta,Z,r,Delta,Iota,R,W,theta)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
opt_beta = zeros(p,1);
beta = ini_beta;  % % 杩唬鍒濆??
% beta = zeros(p,1);
opt_BIC = 1e+10;
n0 = sum(Delta);
lambda0 = log(n0);

for j = 1:length(theta)
    k=1;err=0; tk = 4;
    while k<=1000 && err==0
        %k
        W1 = diag(lambda0*theta(j)*exp(-theta(j)*beta.^2));
        u = eye(p) + 2*W1/tk;

        for i=1:n
            tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
            f2(i,:) = sum(tem1.*Z)/sum(tem1);
        end   

        L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

        beta_tilde = beta - L_prime/tk;
        beta1 = u\beta_tilde;
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end
   
    beta2 = beta.*(abs(beta)>=2*1e-4);

    for i=1:n
        tem2(i,1) = sum((W(:,i).*R(:,i)).*exp(Z*beta2));
    end

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
%                 lasso_lqa()
% ============================================================
function [opt_beta,opt_lambda] = lasso_lqa(n,ini_beta,Z,r,status,Iota,R,W,lambda)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
beta = ini_beta;
opt_BIC = 1e+10;

for i = 1:length(lambda)
    k = 1;err = 0; tk = 4; % beta = ini_beta;
    while k<=1000&&err==0
        %K
        for j=1:n
            tem1 = (W(:,j).*R(:,j)).*exp(Z*beta);
            f2(j,:) = sum(tem1.*Z)/sum(tem1);
        end   
        L_prime = -sum((status.*Iota).*(Z-f2))'/n;
        
        u = beta - L_prime/tk;
        beta1 = max(zeros(p,1),abs(u)-lambda(i)/tk).*sign(u);
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end
    
    beta2 = beta.*(abs(beta)>=2*1e-4);

    for t=1:n
        tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
    end

    ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(i);
        opt_beta = beta2;
    end
end

end

%% ===================================================
%                 Alasso_lqa()
% ============================================================
function [opt_beta,opt_lambda] = Alasso_lqa(n,ini_beta,Z,r,status,Iota,R,W,lambda)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
WW = 1./(abs(ini_beta)+1e-6);
beta = ini_beta;
opt_BIC = 1e+10;
% beta3 = zeros(p,1);

for i = 1:length(lambda)
    k = 1;err = 0;tk = 4; % beta = ini_beta;
    while k<=1000&&err==0
        %       k
        % % === solving by Lasso form ===
        for j=1:n
            tem1 = (W(:,j).*R(:,j)).*exp(Z*beta);
            f2(j,:) = sum(tem1.*Z)/sum(tem1);
        end   
        L_prime = -sum((status.*Iota).*(Z-f2))'/n;

        u = beta - L_prime/tk;
        beta1 = max(zeros(p,1),abs(u)-lambda(i).*WW/tk).*sign(u);      
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k=k+1;
    end
    
    beta2 = beta.*(abs(beta)>=2*1e-4);

    for t=1:n
        tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
    end

    ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(i);
        opt_beta = beta1;
    end
    
end

end



%% ===================================================
%                 scad_lqa()
% ============================================================
function [opt_beta,opt_lambda] = scad_lqa(n,ini_beta,Z,r,status,Iota,R,W,lambda)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
beta = ini_beta;
opt_BIC=1e+10; a = 3.7;

for j=1:length(lambda)
    k=1;err=0; tk = 2; % beta = ini_beta;
    while k<=1000 && err==0
%         k
        W1 = diag( ( lambda(j)*(abs(beta)<=lambda(j)) +...
            ((lambda(j)<abs(beta))&(abs(beta)<a*lambda(j))).*((a*lambda(j)-abs(beta))/(a-1)) )./(abs(beta)+1e-6) );
        u = eye(p) + W1/tk;

        for i=1:n
            tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
            f2(i,:) = sum(tem1.*Z)/sum(tem1);
        end   
        L_prime = -sum((status.*Iota).*(Z-f2))'/n;

        beta_tilde = beta - L_prime/tk;
        beta1 = u\beta_tilde;
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end
    
    beta2 = beta.*(abs(beta)>=2*1e-4);

    for t=1:n
        tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
    end

    ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;
    
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(j);
        opt_beta = beta2;
    end
    
end

end


%% ===================================================
%                 MCP_lqa()
% ============================================================
function [opt_beta,opt_lambda] = MCP_lqa(n,ini_beta,Z,r,status,Iota,R,W,lambda)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
beta = ini_beta;
tv = 2.7;    % according to Cao et al.(2017)
opt_BIC=1e+10;

for i=1:length(lambda)
    k=1;err=0; tk = 4; % beta = ini_beta;
    while k<=1000 && err==0
%        k
        W1=diag(max(0,(tv*lambda(i)-abs(beta))/(tv*lambda(i)) )./(abs(beta)+1e-6));
        u = eye(p) + W1/tk;

        for j=1:n
            tem1 = (W(:,j).*R(:,j)).*exp(Z*beta);
            f2(j,:) = sum(tem1.*Z)/sum(tem1);
        end   

        L_prime = -sum((status.*Iota).*(Z-f2))'/n;

        beta_tilde = beta - L_prime/tk;
        beta1 = u\beta_tilde;
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end
    
    beta2 = beta.*(abs(beta)>=2*1e-4);

    for t=1:n
        tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
    end

    ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(i);
        opt_beta = beta2;
    end
end

end


%% ===================================================
%                 SELO_lqa()
% ============================================================
function [opt_beta,opt_lambda] = SELO_lqa(n,ini_beta,Z,r,status,Iota,R,W,lambda)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
opt_beta = zeros(p,1);
beta = ini_beta ;
% tau = 1.5*0.01;
tau = 0.01;    % according to Dicker et al.(2013)
opt_BIC=1e+10;

for i=1:length(lambda)
    k=1;err=0; tk = 4; % beta = ini_beta;
    while k<=1000 && err==0
        W1=diag((lambda(i)*tau/(2*log(2)))./((abs(beta)+1e-6).*(2*beta.^2+3*tau*abs(beta)+tau^2)));
        u = eye(p) + 2*W1/tk;

        for j=1:n
            tem1 = (W(:,j).*R(:,j)).*exp(Z*beta);
            f2(j,:) = sum(tem1.*Z)/sum(tem1);
        end   

        L_prime = -sum((status.*Iota).*(Z-f2))'/n;

        beta_tilde = beta - L_prime/tk;
        beta1 = u\beta_tilde;
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end

    beta2 = beta.*(abs(beta)>=2*1e-4);

    for t=1:n
        tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
    end

    ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(i);
        opt_beta = beta2;
    end
end

end




%% ===================================================
%                 MIC_lqa()
% ============================================================
function [opt_beta,opt_a] = MIC_lqa(n,ini_beta,Z,r,status,Iota,R,W,a)
[~,p] = size(Z);
f2 = zeros(n,p);
beta = ini_beta;
n0 = sum(status);
lambda0 = log(n0);
% a = 70;
%a = n0;
opt_BIC=1e+10;

for i=1:length(a)
k=1;err=0; tk = 2;
while k<=1000 && err==0
    %   k
    W1 = diag( (4*a(i)*lambda0*exp(2*a(i)*beta.^2))./((exp(2*a(i)*beta.^2)+1).^2) );
    u = eye(p) + 2/tk*W1;
    for j=1:n
        tem1 = (W(:,j).*R(:,j)).*exp(Z*beta);
        f2(j,:) = sum(tem1.*Z)/sum(tem1);
    end  

    L_prime = -sum((status.*Iota).*(Z-f2))'/n;

    beta_tilde = beta - L_prime/tk;
    beta1 = u\beta_tilde;
    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

 beta2 = beta.*(abs(beta)>2*1e-4);

 for t=1:n
     tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
 end

 ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
 sel = beta2~=0;
 BIC = ell+sum(sel)*log(n)/n;
 if BIC<=opt_BIC
     opt_BIC = BIC;
     opt_a = a(i);
     opt_beta = beta2;
 end
end


end


