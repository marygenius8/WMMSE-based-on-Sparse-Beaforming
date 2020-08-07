%%
% 20190209-0210_code_version:
% test the result of this WMMSE_sum_rate function
% clc
% clear
% clear all
% K = 10;
% num_H = 1;
% disp('####### Generate Channel #######');
% generate(K,num_H);
% 
% disp('####### Evaluate Performance #######');
% rate_wmmse_sum=[];
% for loop=1:num_H
%     
%     temp_H = X(:,loop);
%     H=reshape(temp_H,K,K);
% %     Pmax = Pm*ones(K,1);    
%     rate_wmmse = obj_IA_sum_rate(H, Y(:,loop), var_noise);
%     rate_wmmse_sum=[rate_wmmse_sum rate_wmmse];
%     
% end
% % 
% % disp('####### Evaluate Testing Performance #######');
% % testperformance(K,num_H)


%%
% 20190211-0225_code_version:
% the frame of CDF performance dynamic clustering
%%

%%
% 20190306-0314_code_version:
% Algorithm 1 MIMO adaptation and setup communication environment, received signal power adaptation
%%

%%
% 20190316-0319_code_version:
% Check bug - Convergence show!
%%
% 20190407_code_version:
% Check - Why not converged?
% % shrink the problem dimension
% % restore the channel coefficients that dynamic algorithm not converged underlied
%%

% main function to solve joint optimization of SISO Downlink User Association and Power Allocation
clc
clear
clear all
%%%%%%%%%%%%%Problem Simulation Setup%%%%%%%%%%%%%
params.K = 10;         % =# users in each cell
params.L_MB = 1;       % =# Macro base station in each cell
params.L_pB = 1;       % =# Pico cell in each cell
params.M_MB = 2;       % =# antennas at each Macro base station
params.M_pB =2;       % =# antennas at each Pico-cell
params.N_user = 2;     % =# antennas at each user
params.BW = 10;        % bandwidth 10MHz
params.R = 400/sin(pi/3);    % distance between cell 0.8km, radius of each cell 400m
params.center = [0;0];       %Macro base station position (center of cell)
params.Pmax_MB = 19.952;     % the maximum transmit power(Watt) of each Macro base station 
params.Pmax_pB = 1;          % the maximum transmit power(Watt) of each Pico-cell
params.Cmax_MB = 245;     % unit bandwidth capacity
params.Cmax_pB = 70;
% params.var_noise = 10^(-16.9)*10^4;     % 10MHz noise power(Watt)  
params.var_noise = 10^(-16.9)*10^4;     % 10MHz noise power(Watt)  
params.tau = 1e-10;
num_H = 1;            % total number of channel coefficients samples
copy = 1;
% K = 10;
% M = 4;
% N = 2;
% L = 4;
% Pm = 1;
%%%%%%%%%%%%%%%%%generate channel%%%%%%%%%%%%%%%%%
disp('####### Generate channel coefficient #######');
for loop=1:num_H   
    H(:,:,loop)=channel_realization(params);
end
params.H = H;
% params.H = channel_realization(params);

disp('####### Solve WMMSE beamformer #######');
generate_sparse_BF(params,copy,num_H);

disp('####### Load Optimized Data #######');
% load(sprintf('%d_realization%d_user%d_Macro%d_pico%d.mat',copy,num_H, params.K, params.L_MB, params.L_pB));
% % load(sprintf('Test_realization%d_user%d_Macro%d_pico%d.mat',num_H, params.K, params.L_MB, params.L_pB));
% disp('####### Load Analysis Data #######');
% load(sprintf('%d_output_for_analysis_user%d.mat',copy,params.K));

% disp('#### Check Starts ####');
% Converge_init_check(params,num_H);
% disp('#### Check Ends ####');

% disp('####### Evaluate Performance #######');
% rate_wmmse = Z;
% %  CDF of dynamic clustering
% x = sort(10*rate_wmmse(:));
% 
% figure(2);
% xlabel('User Rate (Mbps)');
% ylabel('Cumulative Distribution');
% % plot(x, prob, 'b-o');
% cdfplot(x);
% hold on;
% grid on;
% hold off;
% legend('Algorithm 1 Dynamic Clustering');