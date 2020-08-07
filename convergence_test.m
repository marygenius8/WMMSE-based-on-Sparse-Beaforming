for i=1:10
    flag(i) = Sparse_BF_Clustering_Dynamic_A1_Optimality_test_func(i);
end
converge_num = sum(flag);
converge_rate = sum(flag)/10;
disp('####### Converge Realizations Counts#######');
disp(converge_num);
disp('####### Ratio of Converged Realization/total 10 realizations#######');
disp(converge_rate);