function Converge_init_check(params,num_H)
% copy = 2;
% copy = 3;
copy = 4;
tic
% rng('default');
L = params.L_MB + params.L_pB; LM = params.L_MB * params.M_MB + params.L_pB * params.M_pB;
% X=zeros(K*N,LM,num_H);       % restore channel coefficients matrix under each channel realization 
XX = params.H;
YY = zeros(LM,params.K,num_H);         % restore optimized beamformer matrix under each channel realization
ZZ = zeros(params.K,num_H);            % restore optimized rate vector under each channel realization

% initial beamforming vector  
% bb_int_MB = complex(randn(params.M_MB,params.K),randn(params.M_MB,params.K)); 
% bb_int_MB = bb_int_MB/norm(bb_int_MB,'fro')*sqrt(params.Pmax_MB); 
% bb_int_pB = complex(randn(params.M_pB,params.K),randn(params.M_pB,params.K));
% bb_int_pB = bb_int_pB/norm(bb_int_pB,'fro')*sqrt(params.Pmax_pB); 
cf = 0.8; % change the initial points to check the convergence by manipulating this coefficients
bb_int_MB = complex(ones(params.M_MB,params.K),ones(params.M_MB,params.K)); 
bb_int_MB = cf*bb_int_MB/norm(bb_int_MB,'fro')*sqrt(params.Pmax_MB); 
bb_int_pB = complex(ones(params.M_pB,params.K),ones(params.M_pB,params.K));
bb_int_pB = cf*bb_int_pB/norm(bb_int_pB,'fro')*sqrt(params.Pmax_pB);
bb_int = [];
for l = 1 : L
    if l <= params.L_MB
        bb_int = [bb_int;bb_int_MB];
    else
        bb_int = [bb_int;bb_int_pB];
    end
end
for loop = 1:num_H
    CH = XX(:,:,loop);
%     b_int = 1/sqrt(2)/4*complex(rand(LM,K),rand(LM,K)); % initial beamforming vector     
    [bb_opt,rate_opt,P_iter,SumRate_iter,~, ~] = WMMSE_sparse_BF(params, bb_int, CH,copy);
%     YY(:,:,loop)=b_opt;          % restore optimized beamformer in Y
%     ZZ(:,loop) = rate_opt;       % restore optimized user rate in Z
%     VV(:,:,:,loop) = P_iter;        % restore optimized per BS power per user in final WMMSE iteration:d1=lth BS,d2=kth user,d3=iteration,d4=loop
%     WW(:,loop) = SumRate_iter; % restore sum rate in each WMMSE iteration:d1=sum(rate(1:k)),d2=loop
    if mod(loop,5)==0
        fprintf('.');
    end
end

toc
end