function [train_result_matrix, DistSample_frst_norm_b,DistSample_frst_norm_m] =finddist_splitsets(train_result_matrix,Train_Set,QQ_trans_opt_b,QQ_trans_opt_m,min_val_b_opt,max_val_b_opt,min_val_m_opt,max_val_m_opt,threshold_dist_b,threshold_dist_m)

numRows_tr = size(Train_Set, 2);
resCol_tr_splitb = zeros(numRows_tr, 1);
resCol_tr_splitm = zeros(numRows_tr, 1);
%res_split= zeros(numRows_tr, 1);
train_result_matrix = [train_result_matrix, resCol_tr_splitb, resCol_tr_splitm];%,res_split];
b_col=size((train_result_matrix),2)-1;
m_col=size(train_result_matrix,2);
DistSample_firstsplit_b=[];
DistSample_frst_norm_b=[];
DistSample_firstsplit_m=[];
DistSample_frst_norm_m=[];

%%%%%%%%%%
% There are 2 splitting points available
% 1.  b l1 subspace koptimal distance threshold
% 2.  m's l1 subspace kth rank optimal outlier distance threshold
%
% 1. splitting on b, Q subspace of optrank
%    Finding distance of all records and normalizing it as per training records
DistSample_firstsplit_b=sum(abs(Train_Set - QQ_trans_opt_b* Train_Set),1);
DistSample_frst_norm_b=(DistSample_firstsplit_b-min_val_b_opt)/(max_val_b_opt-min_val_b_opt);

% 2. splitting on m, Q subspace of optrank
%    Finding distance of all records and normalizing it as per training records
DistSample_firstsplit_m=sum(abs(Train_Set - QQ_trans_opt_m* Train_Set),1);
DistSample_frst_norm_m=(DistSample_firstsplit_m-min_val_m_opt)/(max_val_m_opt-min_val_m_opt);

%%%%%Classifying the records as b or n by comparing  their distance with
%%%%%threshold obtained from the algo
for count = 1: size(Train_Set,2)
     
    if DistSample_frst_norm_b(count) < threshold_dist_b
        
        train_result_matrix(count,b_col)=2;
    end
  
    if DistSample_frst_norm_m(count) < threshold_dist_m
        train_result_matrix(count,m_col)=4;
    end
end
% train_result_matrix;
%DistSample_frst_norm_b
% DistSample_frst_norm_m;
end