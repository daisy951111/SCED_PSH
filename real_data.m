format short
clc;clear;close;

%************************** Reading covariates matrix 'Z' *******************************
filename = './real data analysis/GSE5479_Final_processed_data_1.txt';
fileID = fopen(filename);
data1=textscan(fileID,['%*s',repmat('%n',[1,229])],'HeaderLines',1);
fclose(fileID);

filename = 'C:/Users/HUAWEI/Desktop/SH-SEL0/real data analysis/GSE5479_Final_processed_data_2.txt';
fileID = fopen(filename);
data2 = textscan(fileID,['%*s',repmat('%n',[1,175])],'HeaderLines',1);
fclose(fileID);

data1 = cell2mat(data1);
data2 = cell2mat(data2);
Z = [data1,data2]';

%******* Reading Supplemental file data, and changing the 'Sample_name' into 'char' ***********
data = readmatrix('C:/Users/HUAWEI/Desktop/SH-SEL0/real data analysis/Supplemental_file_1.xls');
[~,~,supply]=xlsread('C:/Users/HUAWEI/Desktop/SH-SEL0/real data analysis/Supplemental_Name.xlsx');
for i=1:length(supply)
    name{i,1}=num2str(supply{i,1});
end

%******************** rewritting the data，cancling 'Age' line when the value equals 0 **************************
index=find(isnan(data(:,26))==1);

data(index,:)=[];
name(index,:)=[];
[~,~,temp_grade]=xlsread('C:/Users/HUAWEI/Desktop/SH-SEL0/real data analysis/Grade.xlsx');

for i=1:length(temp_grade)
    if strcmp(temp_grade{i,1},'HIGH')
        grade(i,1)=3;
    elseif strcmp(temp_grade{i,1},'PUNLMP')
        grade(i,1)=2;
    elseif strcmp(temp_grade{i,1},'LOW')
        grade(i,1)=1;
    else
        grade(i,1)=NaN;
    end
end
        
%********************* 找出Dyrskjot用过的数据 ***************************
index1=find(data(:,27)==1|data(:,31)==1);
X=data(index1,20);  %在能用的数据中提取生存时间
V=data(index1,26);  %在能用的数据中提取exposure变量V
Name=name(index1,:);


%************************** Computing 'Delta' ***********************************
temp=data(:,17)+data(:,22);
temp=temp(index1,:);
delta=zeros(length(index1),1);
Delta=delta+(temp~=0);

%************************** Computing 'Iota' ***********************************
I1=find((data(index1,17)==1&data(index1,22)~=1)|(data(index1,17)==2&data(index1,22)~=1)|(data(index1,22)==1));  % had progression or death from bladder cancer
I2=find(data(index1,17)==0&data(index1,22)~=1&data(index1,22)~=0);  % died from other causes.
Iota=zeros(length(index1),1);
Iota(I1,1)=1;
Iota(I2,1)=2;

%****************** 提取基因数据中的样本名称 ***************************
filename = 'C:/Users/HUAWEI/Desktop/SH-SEL0/real data analysis/GSE5479_clinical_information.txt';
fileID = fopen(filename);
sample_name=textscan(fileID,'%s %*[^\n]','HeaderLines',1);
fclose(fileID);
gene=sample_name{1,1};
I=2:2:808;
gene(I,:)=[];

%************************** 匹配数据 ***********************************
I=[];
J=[];
for i=1:length(Name)
    for j=1:length(gene)
        if strcmp(Name{i,1},gene{j,1})
            I=[I;i];
            J=[J;j];
        end
    end
end

Data=[X,Delta,Iota,V];


%************************** Getting all data ***********************************
All_Data=[];%All_Data=[X,Delta,Iota,Z];
for i=1:length(Name)
    All_Data=[All_Data;Data(I(i),:),Z(J(i),:)];
end
X=All_Data(:,1);
Delta=All_Data(:,2);
Iota=All_Data(:,3);
V=All_Data(:,4);
Z=All_Data(:,5:end);


%*************************** Computing 'R' ************************************
n=length(Name);
R=zeros(n,n);
for i=1:n
    for j=1:n
        if X(j,1)>=X(i,1)||(X(j,1)<=X(i,1)&&Delta(j,1)==1&&Iota(j,1)~=1)
            R(i,j)=1;
        end
    end
end


%************************** Computing 'W' ***********************************
T=X;
C=zeros(length(Name),1);
C(Delta==0,1)=X(Delta==0,1);
%W=Weight(X,T,C,length(Name));
W=Weight(X,T,C,Delta,length(Name));
% W = Weight(X,T,C,Delta,n);



%% === 'theta' & 'delta' parameters ===
% theta2 = (0.20:-0.005:0.055)*1.25*10^4;  % 1000左右
% theta2 = 2000:100:25000;

theta = (0.90:-0.0095:0.60)*1.25*10^4;
delta = (0.075:-0.0009:0.005);

delta1 = (0.45:-0.05:0.35);

n = length(Name); % sample size
d = 52; % screening size

% Z = normalize(Z);  % 标准化
% Z = zscore(Z); 

%% Method 1 (Selecting significant variables bassed on Tian et al.(2024)'s screening results)
% %  sceening set by Tian et al.(2024)
IN = [758 251 370 709 1339 424 240 767 369 997 765 783 823 368 1357 785 768 1142 820 711 774 979 1354 766 964 1379 347 1191 963 ...
      713 79 1084 698 1355 759 1019 49 1110 980 1091 961 438 1310 1312 185 577 184 1038 710 977 1356 334];

% % Variable selection
initial_beta1 = ini_beta2(Z(:,IN),1e-3,Delta,Iota,R,W);  % % initial beta
[hlselo_beta,opt_theta,opt_delta] = hlselo_lqa(n,initial_beta1,Z(:,IN),1e-4,Delta,Iota,R,W,theta,0.05);

index2 = find(hlselo_beta~=0);

% IN(index2)
% [370 709 767 369 713] %  by Tian, Liu and Wang (2022) VCR-IGHT

cov_beta = cov_PSH(index2,Delta,Iota,Z(:,IN),W,R,hlselo_beta);  % covariance of coefficient

PSH_std_beta  = sqrt(cov_beta);
z_stat1 = hlselo_beta(index2)./PSH_std_beta;
p_value1 = (1-normcdf(abs(z_stat1)))*2;   %% p_value

[hlselo_beta(index2) PSH_std_beta p_value1]



%% Method 2
% % Feature screening
initial_beta2 = ini_beta1(Z,1e-3,Delta,Iota,R,W);  % % initial beta

opt_beta = NonMargScr_PSH(n,initial_beta2,Z,1e-3,Delta,Iota,R,W,d); %  % feature screening

index3 = find(opt_beta~=0);

initial_beta3 = ini_beta2(Z(:,index3),1e-3,Delta,Iota,R,W);  % % initial beta  
[hlselo_beta2,opt_theta2,opt_delta2] = hlselo_lqa(n,initial_beta3,Z(:,index3),1e-4,Delta,Iota,R,W,theta,delta);  % % HLSELO

index4 = find(hlselo_beta2~=0);

cov_beta2 = cov_PSH(index4,Delta,Iota,Z(:,index3),W,R,hlselo_beta2);  % covariance of coefficient

PSH_std_beta2  = sqrt(cov_beta2);
z_stat2 = hlselo_beta2(index4)./PSH_std_beta2;
p_value2 = (1-normcdf(abs(z_stat2)))*2;   %% p_value

[hlselo_beta2(index4) PSH_std_beta2 p_value2]


%% Method 3
block_size = 20;
initial_beta_cd = ini_beta_cd(Z,1e-3,Delta,Iota,R,W,theta,block_size);  % % initial beta - CD algorithm  
block_size
[hlselo_beta_cd,opt_theta3,opt_delta3] = hlselo_CD(n,initial_beta_cd,Z,1e-4,Delta,Iota,R,W,theta,delta1,block_size);

% [hlselo_beta_cd,opt_theta3,opt_delta3] = hlselo_CD(n,initial_beta_cd,Z,1e-4,Delta,Iota,R,W,theta,0.055,block_size);
index5 = find(hlselo_beta_cd~=0);

Z3 = Z(:,index5);
beta3 = hlselo_beta_cd(index5);

cov_beta3 = cov_PSH(1:length(index5),Delta,Iota,Z3,W,R,beta3);

PSH_std_beta3 = sqrt(cov_beta3);
z_stat3 = beta3./PSH_std_beta3;
p_value3 = (1-normcdf(abs(z_stat3)))*2;

[beta3 PSH_std_beta3 p_value3]


% % %************************** 读取所有Prob名称 ***********************************
supply_name = readcell('C:/Users/HUAWEI/Desktop/SH-SEL0/real data analysis/GSE5479_Final_processed_data_1.txt');
name_prob = string(supply_name(2:1382,1));

% %************************** 读取第二阶段的重要Prob名称 *************************
selected_prob_name1 = name_prob(index2) % % Method 1

selected_prob_name2 = name_prob(index4) % % Method 2

selected_prob_name3 = name_prob(index5) % % Method 3




%% ============================================================
%                  SUBFUNCTIONS
% ===============%==%==%==%==%==%==%===========================


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
% [f,x] = ecdf(C,'Censoring',Delta); %删失时间的KM估计 
% G = 1-f;
% ecdf(C,'censoring',Delta,'function','survivor'); %删失时间的生存函数 

[n,~] = size(X);
cens = 1-Delta;
[~,id] = sort(X);
cens = cens(id);  % 删失指示变量

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
%                 ini_beta1()
% ============================================================
function initial_beta = ini_beta1(Z,r,Delta,Iota,R,W)
[n,p] = size(Z);
beta = zeros(p,1);
f2 = zeros(n,p);

k = 1;err = 0; tk = 280;
while k<=1000&&err==0
    %k
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

initial_beta = beta;

end



%% ===================================================
%                 ini_beta2()
% ============================================================
function initial_beta = ini_beta2(Z,r,Delta,Iota,R,W)
[n,p] = size(Z);
beta = zeros(p,1);
f2 = zeros(n,p);
k = 1;err = 0; tk = 120; % tk = 240;
while k<=1000&&err==0
    %k
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
initial_beta = beta;

end




%% ===================================================
%                 ini_beta_cd()
% ============================================================
function initial_beta = ini_beta_cd(Z,r,Delta,Iota,R,W,theta,block_size)
[n,p] = size(Z);
beta = zeros(p,1);
% f2 = zeros(n,p);

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


k = 1;err = 0; tk = 20;
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

% for j = 1:length(theta)
% 
% k=1;err=0; tk = 16;
% while k<=1000 && err==0
%     %k
%     W1 = lambda0*theta(j)*diag(1./(1+theta(j)*beta.^2).^2);
%     u = eye(p) + 2*W1/tk;
% 
%     % for i=1:n
%     %     tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
%     %     %tem1 = R(:,i).*exp(Z*beta);
%     %     f2(i,:) = sum(tem1.*Z)/sum(tem1);
%     % end
% 
%     % revised on 17/03/2025
%     tem1 = (W.*R).*exp(Z*beta);
%     sum_tem1 = sum( tem1, 1 );
%     f2 = (Z'*tem1)./sum_tem1;
%     f2 = f2';
% 
%     L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;
% 
%     beta_tilde = beta - L_prime/tk;
%     beta1 = u\beta_tilde;
%     w = beta1-beta;
%     err = norm(w,2)^2 <= r*norm(beta,2)^2;
%     beta = beta1;
%     k = k+1;
% end
% 
% beta2 = beta.*(abs(beta)>=2*1e-4);
% 
% % for i=1:n
% %     tem2(i,1) = sum((W(:,i).*R(:,i)).*exp(Z*beta2));
% % end
% 
% % revised on 17/03/2025
% tem2 = sum( (W.*R).*exp(Z*beta),1 );
% 
% 
% ell = -sum((Delta.*Iota).*(Z*beta2-log(tem2)))/n;
% sel = beta2~=0;
% BIC = ell+sum(sel)*log(n)/n;
% 
% if BIC<=opt_BIC
%     opt_BIC = BIC;
%     opt_beta = beta2;
% end
% 
% end
% 
% initial_beta = opt_beta;



end




%% ===================================================
%                 hlselo_lqa()
% ============================================================
function [opt_beta,opt_theta,opt_delta] = hlselo_lqa(n,ini_beta,Z,r,Delta,Iota,R,W,theta,delta)
[~,p] = size(Z);
% f2 = zeros(n,p);
% tem2 = zeros(n,1);
opt_beta = zeros(p,1);
beta1 = zeros(p,1);
% beta = ini_beta;  % % initial value
% beta = zeros(p,1);
opt_BIC = 1e+10;
n0 = sum(Delta);
lambda0 = log(n0);

for h = 1:length(delta)
    for j = 1:length(theta)
        beta = ini_beta;
        k=1;err=0; tk = 300; % beta = ini_beta;
        while k<=1000 && err==0
            %k

            tem1 = (W.*R).*exp(Z*beta);
            sum_tem1 = sum( tem1, 1 );
            f2 = (Z'*tem1)./sum_tem1;
            f2 = f2';

            L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

            beta_tilde = beta - L_prime / tk;
            index = find(abs(beta_tilde) > lambda0*delta(h));
            index1 = setdiff(1:p,index);

            beta1(index1) = zeros(length(index1),1);
            W1 = lambda0 * (2 * abs(beta(index)) * theta(j) + delta(h)) .* ...
                diag(1 ./ (1 + delta(h)*abs(beta(index)) + theta(j)*beta(index).^2).^2 + 1e-6);
            u = eye(length(index)) + 2 * W1 / tk;
            beta1(index) = u \ beta_tilde(index);

            w = beta1-beta;
            err = norm(w,2)^2 <= r*norm(beta,2)^2;
            beta = beta1;
            k = k+1;
        end

        beta2 = beta1;

        tem2 = sum( (W.*R).*exp(Z*beta2),1 );
        ell = -sum((Delta.*Iota).*(Z*beta2-log(tem2)))/n;
        sel = beta2~=0;
        % sel = beta2 > 1e-6;
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
%                 hlselo_CD()
% ============================================================
function [opt_beta,opt_theta,opt_delta] = hlselo_CD(n,ini_beta,Z,r,status,Iota,R,W,theta,delta,block_size)
[~,p] = size(Z);
opt_beta = zeros(p,1);
% beta = ini_beta;  % % Initial iteration value
beta1 = zeros(p,1);
opt_BIC = 1e+10;
n0 = sum(status);
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

% blocks

% beta_new1 = zeros(block_size,1);

for h = 1:length(delta)
   
    for j = 1:length(theta)
        k=1; err=0; tk = 50; beta = ini_beta;
        while k<=1000 && err==0
            %k

            tem1 = (W.*R).*exp(Z*beta);
            sum_tem1 = sum( tem1, 1 );
            f2 = (Z'*tem1)./sum_tem1;
            f2 = f2';

            L_prime = -sum((status.*Iota).*(Z-f2))'/n;

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

                        block_len = length(current_block);  % 动态获取当前块长度
                        beta_block = beta(current_block);

                        % 动态创建 beta_new1，大小匹配当前块
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

                    err = norm(beta - beta1,2)^2 / (norm(beta,2)^2 + 1e-8);

                    if err < 1e-6
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
        ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
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
%                 lselo_LQA()
% ============================================================
function [opt_beta,opt_theta] = lselo_LQA(n,ini_beta,Z,r,Delta,Iota,R,W,theta)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
opt_beta = zeros(p,1);
beta = ini_beta;  % % 迭代初值
% beta = zeros(p,1);
opt_BIC = 1e+10;
n0 = sum(Delta);
lambda0 = log(n0);

%theta = 1000;

for j = 1:length(theta)
    k=1;err=0; tk = 500; %beta = ini_beta; tk = 500;
    while k<=1000 && err==0
        k
        W1 = lambda0*theta(j)*diag(1./(1+theta(j)*beta.^2).^2);
        u = eye(p) + 2*W1/tk;

        for i=1:n
            tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
            %tem1 = R(:,i).*exp(Z*beta);
            f2(i,:) = sum(tem1.*Z)/sum(tem1);
        end   

        L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

        beta_tilde = beta - L_prime/tk;
%         beta_tilde = beta - L_prime/(n*tk);
        beta1 = u\beta_tilde;
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end
   
    beta2 = beta.*(abs(beta)>=3*1e-4);

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
%                 NonMargScr_PSH()
% ============================================================
function opt_beta = NonMargScr_PSH(n,ini_beta,Z,r,Delta,Iota,R,W,d)
[~,p] = size(Z);
beta = ini_beta;
% beta = zeros(p,1);

k=1;err=0; tk = 104;
while k<=1000 && err==0
    %k
     for i=1:n
            tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
            %tem1 = R(:,i).*exp(Z*beta);
            f2(i,:) = sum(tem1.*Z)/sum(tem1);
     end   

    L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

    beta_tilde = beta - L_prime/tk;

    b2 = sort(abs(beta_tilde),'descend');
    beta1 = beta_tilde.*(abs(beta_tilde) >= b2(d));

    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;

end

opt_beta = beta1;

end



%% ===================================================
%                 cov_PSH()
% ============================================================
%function cov_beta = cov_PSH(Z,r,Delta,Iota,R,W);

function cov_beta = cov_PSH(index,Delta,Iota, Z, W, R ,beta)
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
  f1 = f1 + ( DI(i)/tt(i) )*Cel3{1,i};   % % 二阶偏导的第一项
  f2 = f2 + ( DI(i)/tt(i)^2 )*Cel4{1,i}; % % 二阶偏导的第二项
end


% f1 = cellfun(@sum, Cel1, 'UniformOutput', false);   % % 二阶偏导的第一项
% f2 = cellfun(@sum, Cel2, 'UniformOutput', false);   % % 二阶偏导的第二项


L_primeprime = (-f1+f2);

cov_beta = diag( inv( L_primeprime(index,index) ) );  


end