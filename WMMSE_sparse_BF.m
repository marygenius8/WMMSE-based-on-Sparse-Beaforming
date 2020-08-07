%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     WMMSE approach to generate optimized beamformer matrix
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b_opt, r_opt, P_iter,SumRate_iter,sr_opt, totIters] = WMMSE_sparse_BF(params,b_int, H,copy)
% initialize every variables
K = params.K; L = params.L_MB + params.L_pB; M_MB = params.M_MB; M_pB = params.M_pB;
N = params.N_user; LM = params.L_MB * M_MB + params.L_pB * M_pB; var_noise = params.var_noise;
tau = params.tau;          % small constant regularization factor used for beta update
rnew = 0;             % sum rate
b = b_int;            % (w_k) square root of the transmit power as initial power 
u = complex(zeros(N,K));    % use u to restore receiver u_k, receiver beformer or amplifier coefficient?
rho = zeros(K, 1);    % rho is iterative updating weight
beta = zeros(K, L);   % beta is iterative updating weight
a = zeros(K, 1);      % alpha, weights of the WSR rates
a(:) = 1;             %  
Pmax = [params.Pmax_MB*ones(params.L_MB,1);params.Pmax_pB*ones(params.L_pB,1)];
Cmax = [params.Cmax_MB*ones(params.L_MB,1);params.Cmax_pB*ones(params.L_pB,1)];
Rate = zeros(K,1);
interf = zeros(N,N,K); % interference
interf_noise = zeros(N,N,K); % sum of interference and noise
SINR = zeros(K, 1);
I = eye(N);           % identity matrix used for interference
mv01 = zeros(N,1);     % mv and m1v is the intermediate variables used for interference calculate
mv02 = zeros(N,N);      
mv1 = zeros(LM,N);
mv21 = zeros(N,1);
mv22 = zeros(LM,N);
% B_iter = [] ;% record the beamforming vector in each WMMSE iteration

%% initialize the u_0,rho_0,rnew(sum rate)
for i=1 : K
    idx_lpB = params.L_MB * params.M_MB;
    for j = 1 : K                      % update SINR        
        if j ~= i
            interf(:,:,i) = interf(:,:,i) + H(N*(i-1)+1:N*i,:)*b(:,j)*(H(N*(i-1)+1:N*i,:)*b(:,j))';
        end
    end
    sig = H(N*(i-1)+1:N*i,:)*b(:,i);
    interf_noise(:,:,i) = interf(:,:,i) + var_noise*I;    % interference plus noise
%     SINR(i) =  sig'/((interf_noise(:,:,i)+conj(transpose(interf_noise(:,:,i))))/2)*sig; % current SINR
    SINR(i) =  sig'/interf_noise(:,:,i)*sig; % current SINR
    Rate(i) =  10*log2(1+real(SINR(i)));          % Rate (k)
%     SINR(i) = conj(H(i,:)*b(:,i))*(H(i,:)*b(:,i))/interf_noise(i); % current SINR
%     SINR(i) = (abs(H(i,:)*b(:,i)))^2/interf_noise(i); % current SINR
%%     initialize backhual weights beta(j,i)%%
    for l = 1 : L
        if l<=params.L_MB
            beta(i,l) = 1/(norm(b(M_MB*(l-1)+1:M_MB*l,i),'fro')^2+tau);  % Since the outside iteration is the user column, 
%             idx_lpB = l*params.M_MB;
        else
            beta(i,l) = 1/(norm(b(idx_lpB+(1:M_pB),i),'fro')^2+tau);  % Since the outside iteration is the user column, 
            idx_lpB = idx_lpB + params.M_pB;
        end
    end                                   % then the inner iteration is the RRH row when the user k (column k) is fixed

    sig_sq = sig*sig';
    u(:,i) = inv(interf_noise(:,:,i) + sig_sq)*sig;  % update receiver u(i)
%     u(:,i) = inv(((interf_noise(:,:,i) + sig_sq)+conj(transpose(interf_noise(:,:,i) + sig_sq)))/2)*sig;  % update receiver u(i)
    rho(i) = 1/(1+real(u(:,i)'*(interf_noise(:,:,i) + sig_sq)*u(:,i))-2*real(u(:,i)'*sig));  
end

%% iteratively to calculate w_k and update R_k and Beta_k, u(k) and rho(k) optimal until convergence
SumRate = [sum(Rate)];

% Iteration
iter = 0;
while(1)
    iter = iter+1;
    rold = rnew;
%%   calculate beamformer vector b(i) based on QCQP CVX
%     disp('Compute beamformer through QCQP...')
%%   update beamformer b in (LM,K) dimension   %%
    cvx_begin quiet
        obj = 0;
        variable x(LM,K) complex
        for k = 1:K
            obj_sub = zeros(LM,LM);
            for j = 1:K
                mv1 = H(N*(j-1)+1:N*j,:)'*u(:,j);
                obj_sub = obj_sub + a(j) * rho(j) * (mv1 * mv1');
            end
            obj = obj + quad_form(x(:,k), obj_sub ) - 2 * a(k) * rho(k) * real(u(:,k)' * H(N*(k-1)+1:N*k,:) * x(:,k));
%             obj = obj + x(:,k)'*obj_sub* x(:,k)- 2 * a(k) * rho(k) * real(u(:,k)' * H(N*(k-1)+1:N*k,:) * x(:,k));
%             obj = obj + quad_form( x(:,k), obj_sub);    
        end
        minimize(obj)
        subject to
            for l = 1 : L                   %%%Transmit Power Constraints
                if l <= params.L_MB
                    square_pos(norm(x(M_MB*(l-1)+1:M_MB*l,:),'fro')) <= Pmax(l);
                    idx_lpB = l*params.M_MB;
                else
                    square_pos(norm(x(idx_lpB+(1:M_pB),:),'fro')) <= Pmax(l);
                    idx_lpB = idx_lpB + params.M_pB;
                end
            end
            
            for l = 1 : L                   %%% Back hual Constraints
                if l <= params.L_MB
                    square_pos(norm(x(M_MB*(l-1)+1:M_MB*l,:)*sqrt(beta(:,l).*Rate),'fro')) <= Cmax(l);
                    idx_lpB = l*params.M_MB;
                else
                    square_pos(norm(x(idx_lpB+(1:M_pB),:)*sqrt(beta(:,l).*Rate),'fro')) <= Cmax(l);
                    idx_lpB = idx_lpB + params.M_pB;
                end
            end
    cvx_end

    b = x;
%% update variables for iteration to zeros%%
    % update R_i and Beta_i, u(i) and rho(i)
    rnew = 0;
    rate = [ ];  % optimized each user rate
    interf = zeros(N,N,K); % interference
    interf_noise = zeros(N,N,K); % interference + noise
%%  update parameters for each user(i)
    for i = 1 : K
%         interf(i) = (norm(H*b(:,i),2))^2 - (abs(H(i,:)*b(:,i)))^2
%%   update iterference(i)   %% 
        for j = 1 : K                      
            if j~=i
                interf(:,:,i) = interf(:,:,i) + H(N*(i-1)+1:N*i,:)*b(:,j)*(H(N*(i-1)+1:N*i,:)*b(:,j))';
            end
        end
        sig = H(N*(i-1)+1:N*i,:)*b(:,i);
%%   update SINR(i)   %%
        interf_noise(:,:,i) = interf(:,:,i) + var_noise*I;    
%         SINR(i) = mv21'/interf_noise(:,:,i)*mv21; 
%         SINR(i) = sig'/interf_noise(:,:,i)*sig;
        SINR(i) = sig'/((interf_noise(:,:,i)+conj(transpose(interf_noise(:,:,i))))/2)*sig;
%%   update Rate(i)   %% 
        Rate(i) =  10*log2(1+real(SINR(i)));         % update the Rate (k)
        % there needs the weights for each users' rate to get weighted sum rate to stop the iteration
        rnew = rnew + a(i)*Rate(i);             % cumulate rnew and update rnew when rate of all user(i)(user_k) updated completely
%%   update backhual weights beta(j,i)   %%
        for l = 1 : L                      
            if l<=params.L_MB
               beta(i,l) = 1/(norm(b(M_MB*(l-1)+1:M_MB*l,i),'fro')^2+tau);  % Since the outside iteration is the user column, \
               idx_lpB = l*params.M_MB;
            else
                beta(i,l) = 1/(norm(b(idx_lpB+(1:M_pB),i),'fro')^2+tau);  % Since the outside iteration is the user column, 
                idx_lpB = idx_lpB + params.M_pB;
            end
        end                                   % then the inner iteration is the RRH row when the user k (column k) is fixed
%%   update receiver u(i)   %%
        sig_sq = sig*sig';% update receiver u(i)
        u(:,i) = inv(interf_noise(:,:,i) + sig_sq)*sig;
%         u(:,i) = inv(((interf_noise(:,:,i) + sig_sq)+conj(transpose(interf_noise(:,:,i) + sig_sq)))/2)*sig;
%%   update receiver rho(i) and rate vector   %%
        rho(i) = 1/(1+real(u(:,i)'*(interf_noise(:,:,i) + sig_sq)*u(:,i))-2*real(u(:,i)'*sig)); % update rho(i)
        rate = [rate;Rate(i)];             % rate of each user
    end
    fprintf('Current iteration is: %d:\n',iter);
    B_iter(:,:,iter) = b;
    SumRate = [SumRate; sum(Rate)]; % record sum capacity evolution  
%%     stop condition
%     if abs(rnew-rold) <= 0.5*1e-2 || iter>=200
    if iter>=1000
%     if abs(rnew-rold) <= 1e-3 || iter>=200
        break;
    end
%     fprint('Each User Rate in each iteration:\n');
%     disp(Rate);
end
B(:,:,[2:iter+1]) = B_iter;
B(:,:,1) = b_int;
%% pick user k's corresponding BS and pico cell TX power  %%
for k = 1:params.K
    idx_lpB = params.L_MB * params.M_MB;
    for l =1:L
        if l <= params.L_MB
            P_iter(l,k,:) = sum(abs(B(params.M_MB*(l-1)+1:params.M_MB*l,k,:)).^2);  % each BS TX power in all iterations
%             idx_lpB = params.M_MB*l;
        else
            P_iter(l,k,:) = sum(abs(B(idx_lpB+(1:params.M_pB),k,:)).^2);
            idx_lpB = idx_lpB+params.M_pB;    
        end
    end
end
% % change the 3D power (D1=lth_BS, D2=kth_user, D3=iteration) to 2D
% matrix,concatinate the iteration dimension
% for ii = 1:iter+1
%     P2D = [P2D reshape(P_iter(:,:,ii),L,params.K) PNaN];
% end
% % capsule BF, Power_iter, rate(1:k) in a 2D matrix output
% outputNaN = nan(1,length(B2D));
% output = [B2D;outputNaN;real(P2D);outputNaN;real(Rate1D)];
outputNaN = nan(1,params.K);
output = [b;outputNaN;SINR';outputNaN;P_iter(:,:,end)];
save(sprintf('%d_output_for_analysis_user%d.mat',copy,params.K),'output');
fprintf('Analysis Data Generate Done! \n');

% draw Fig of the power evolution for each user
% for k = 1:params.K
%     figure(k);
% %     legend('BS1','BS2','BS3','BS4');
%     Pow = reshape(P_iter(:,k,:),L,iter+1)*1e-4;
% %     plot(0:iter,10*log10(Pow),'-o','Displayname',{'BS1','BS2','BS3','BS4'});
%     plot(0:iter,10*log10(Pow),'-o');
%     legend({'BS1','BS2','BS3','BS4'});
%     xlabel(sprintf('k=%d iteration',k));
%     ylabel('Power(dBm/Hz)');
%     title(sprintf('%dMacro-%dpico-%duser',params.L_MB,params.L_pB,params.K));
% end
pick_idx_u = 2;
figure(copy);
Pow = (reshape(P_iter(:,pick_idx_u,:),L,iter+1))*1e-4;    % mw per HZ
plot(0:iter,10*log10(Pow),'-o');
legend({'BS1','BS2','BS3','BS4'});
xlabel(sprintf('k=%d iteration',pick_idx_u));
ylabel('Power(dBm/HZ)');
title(sprintf('%dMacro-%dpico-%duser',params.L_MB,params.L_pB,params.K));

totIters = iter;
sr_opt = rnew;
b_opt = b;
r_opt = rate;
SumRate_iter = SumRate;
fprintf('Total iteration is %d \n',totIters);
fprintf('Each User Rate after total %d iteration:\n',totIters);
disp(rate);

end
