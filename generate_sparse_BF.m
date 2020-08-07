%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate_sparse_BF
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function generate_sparse_BF(params,copy,num_H)

% function generate(prams,num_H,copy)
tic
% rng('default');
L = params.L_MB + params.L_pB; LM = params.L_MB * params.M_MB + params.L_pB * params.M_pB;
% X=zeros(K*N,LM,num_H);       % restore channel coefficients matrix under each channel realization 
X = params.H;
Y = zeros(LM,params.K,num_H);         % restore optimized beamformer matrix under each channel realization
Z = zeros(params.K,num_H);            % restore optimized rate vector under each channel realization

% initial beamforming vector  
% b_int_MB = complex(randn(params.M_MB,params.K),randn(params.M_MB,params.K)); 
% b_int_MB = b_int_MB/norm(b_int_MB,'fro')*sqrt(params.Pmax_MB); 
% b_int_pB = complex(randn(params.M_pB,params.K),randn(params.M_pB,params.K));
% b_int_pB = b_int_pB/norm(b_int_pB,'fro')*sqrt(params.Pmax_pB); 

b_int_MB = complex(ones(params.M_MB,params.K),ones(params.M_MB,params.K)); 
b_int_MB = b_int_MB/norm(b_int_MB,'fro')*sqrt(params.Pmax_MB); 
b_int_pB = complex(ones(params.M_pB,params.K),ones(params.M_pB,params.K));
b_int_pB = b_int_pB/norm(b_int_pB,'fro')*sqrt(params.Pmax_pB);
b_int = [];
for l = 1 : L
    if l <= params.L_MB
        b_int = [b_int;b_int_MB];
    else
        b_int = [b_int;b_int_pB];
    end
end
for loop = 1:num_H
    CH = X(:,:,loop);
%     b_int = 1/sqrt(2)/4*complex(rand(LM,K),rand(LM,K)); % initial beamforming vector     
    [b_opt,rate_opt,P_iter,SumRate_iter,~, ~] = WMMSE_sparse_BF(params, b_int, CH,copy);
    Y(:,:,loop)=b_opt;          % restore optimized beamformer in Y
    Z(:,loop) = rate_opt;       % restore optimized user rate in Z
    V(:,:,:,loop) = P_iter;        % restore optimized per BS power per user in final WMMSE iteration:d1=lth BS,d2=kth user,d3=iteration,d4=loop
    W(:,loop) = SumRate_iter; % restore sum rate in each WMMSE iteration:d1=sum(rate(1:k)),d2=loop
    if mod(loop,5)==0
        fprintf('.');
    end
% save(sprintf('Gussian_%d_realization%d_user%d_Macro%d_pico%d.mat',copy,num_H, params.K, params.L_MB, params.L_pB),'X','Y','Z');
save(sprintf('%d_realization%d_user%d_Macro%d_pico%d.mat',copy,num_H, params.K, params.L_MB, params.L_pB),'X','Y','Z','V','W');
fprintf('Generate Done! \n');
toc
end
