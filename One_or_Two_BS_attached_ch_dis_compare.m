clc
clear
clear all
for copy =1:10
    Sparse_BF_Clustering_Dynamic_A1_Optimality_test_func(copy);
end
% Sparse_BF_Clustering_Dynamic_A1_Optimality_test_func(200);
% figure(200);
% % ColorSet = varycolor(4);  
% % set(gcf, 'Colormap', ColorSet );
% 
idx = [3,7,8]; % ith iteration
for i=1:length(idx)
    load(sprintf('test%d_realization1_user5_Macro1_pico1.mat',idx(i)));
    CH_gain = abs(X);
    for k =1:5
        CH_user(k,:) = sum(CH_gain((k-1)*2+(1:2),:));
    end
    x_ch(i,:) = sort(CH_user(:));
% %     x_ch(i,:) = sort(CH_gain(:));
    cdfplot(x_ch(i,:));
    hold all
end
% 
% % cdfplot(x_ch);
% 
ylim([0 1.5])
xticks(0:0.005:0.1)
xlabel('Channel coefficiets of each user');
ylabel('Cumulative Distribution');
legend({'realization1','realization2','realization3','realization4'})
title('Channel CDF comparation between different realizaiton assignment')
grid on;
hold off;
% 
% figure(201);
% for i=1:4
%     load(sprintf('test%d_output_for_analysis_user5',i));
%     Pow_aloc = abs(output(6:7,:));
%     x_pow(i,:) = sort(Pow_aloc(:));
%     cdfplot(x_pow(i,:));
%     hold all
% end
% 
% ylim([0 1.5])
% % xticks(0:0.005:0.1)
% xlabel('Allocated power of each user');
% ylabel('Cumulative Distribution');
% legend({'realization1','realization2','realization3','realization4'})
% title('Power CDF comparation between different realizaiton assignment')
% grid on;
% hold off;
% 
