clc;clear;close all;
%% 加载数据
theta_true = load('8k512s.mat').air512s;
speech = load('speech.mat').speech; 
para_num = size(theta_true,1);
region_index = 2000;
%% 系统建模
N = 48000;
sound1 = speech(1:N,1);
sound1 = sound1./var(sound1);
sound2 = filter(theta_true,1,sound1);
sound2_two = filter(-theta_true,1,sound1);
sound2 = [sound2;sound2_two];
mu = 0; sigma = 0.1; p = 3; % 广义高斯分布噪声参数
%% 辨识算法+CRLB
for iter = 1:1
    sound2_new = GGD_Model(sound2,2*N,mu,sigma,p);
    input = sound1; theta = zeros(para_num,1);
    % NLMS算法
    d = sound2_new(1:N); alpha = 0.5; delta = 20; [theta_est_NLMS,MSD_NLMS] = algo_NLMS(theta,theta_true,para_num,alpha,delta,input,d); % NLMS算法迭代
    d = sound2_new(N+1:2*N); [theta_est_NLMS_2,MSD_NLMS_2] = algo_NLMS(theta_est_NLMS,-theta_true,para_num,alpha,delta,input,d);
    MSD_NLMS_final = [MSD_NLMS;MSD_NLMS_2(1023:N)]; 
    save(['.\CRLB_Algorithm_Sparse\MSD_NLMS\' num2str(iter) '.mat'],'MSD_NLMS_final');
    % RZA-NLMS算法
    d = sound2_new(1:N); rou_RZA = 1*10^-6; esilon_RZA = 15; alpha = 0.5; delta = 20; [theta_est_RZA_NLMS,MSD_RZA_NLMS] = algo_RZA_NLMS(theta,theta_true,para_num,rou_RZA,esilon_RZA,alpha,delta,input,d); % RZA_NLMS算法迭代
    d = sound2_new(N+1:2*N); [theta_est_RZA_NLMS_2,MSD_RZA_NLMS_2] = algo_RZA_NLMS(theta_est_RZA_NLMS,-theta_true,para_num,rou_RZA,esilon_RZA,alpha,delta,input,d);
    MSD_RZA_NLMS_final = [MSD_RZA_NLMS;MSD_RZA_NLMS_2(1023:N)]; save(['.\CRLB_Algorithm_Sparse\MSD_RZA_NLMS_iter\' num2str(iter) '.mat'],'MSD_RZA_NLMS_final');
    % IPNLMS算法
    d = sound2_new(1:N); mu_IPNLMS = 0.1; alpha_IPNLMS = 0; delta_IPNLMS = 0.001; [theta_est_IPNLMS,MSD_IPNLMS] = algo_IPNLMS(theta,theta_true,para_num,mu_IPNLMS,alpha_IPNLMS,delta_IPNLMS,input,d); % IPNLMS算法迭代
    d = sound2_new(N+1:2*N); [theta_est_IPNLMS_2,MSD_IPNLMS_2] = algo_IPNLMS(theta_est_IPNLMS,-theta_true,para_num,mu_IPNLMS,alpha_IPNLMS,delta_IPNLMS,input,d);
    MSD_IPNLMS_final = [MSD_IPNLMS;MSD_IPNLMS_2(1023:N)]; save(['.\CRLB_Algorithm_Sparse\MSD_IPNLMS\' num2str(iter) '.mat'],'MSD_IPNLMS_final');
    % RLS算法
    d = sound2_new(1:N); lambda_RLS = 0.998; [theta_est_RLS,MSD_RLS] = algo_RLS(theta,theta_true,para_num,lambda_RLS,input,d); % RLS算法迭代
    d = sound2_new(N+1:2*N); [theta_est_RLS_2,MSD_RLS_2] = algo_RLS(theta_est_RLS,-theta_true,para_num,lambda_RLS,input,d); 
    MSD_RLS_final = [MSD_RLS;MSD_RLS_2(1023:N)]; save(['.\CRLB_Algorithm_Sparse\MSD_RLS\' num2str(iter) '.mat'],'MSD_RLS_final');
    % RGM-RLS算法
    d = sound2_new(1:N); lambda_RGM = 0.998; sigma_RGM = 0.3; [theta_est_RGM_RLS,MSD_RGM_RLS] = algo_RGM_RLS(theta,theta_true,para_num,lambda_RGM,sigma_RGM,input,d); % GM_RLS算法迭代
    d = sound2_new(N+1:2*N);  [theta_est_RGM_RLS_2,MSD_RGM_RLS_2] = algo_RGM_RLS(theta_est_RGM_RLS,-theta_true,para_num,lambda_RGM,sigma_RGM,input,d); 
    MSD_RGM_RLS_final = [MSD_RGM_RLS;MSD_RGM_RLS_2(1023:N)]; save(['.\CRLB_Algorithm_Sparse\MSD_RGM_RLS\' num2str(iter) '.mat'],'MSD_RGM_RLS_final');
    % 计算CRLB
    CRLB_history = compute_CRLB(input,para_num,sigma,p); 
end 
