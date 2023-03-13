function [rdos_score,data,data_normal]=outlier_det(data)
%%========================================
% Outlier detection using rdos_outlier
% Input
%-------------
% data   input data Nx3 array 
%
% Output
%------------
% rdos_score               outlierness
% data                     normalized data
% data_normal              record the mean and std of the normal
%%==========================================

%% Normalization of the input data---> zero mean, unit std 
[data,data_normal]=data_normalize_input(data);

[num_data,~]=size(data);



%% the highest k values to examine
K = min(100, num_data-1);


%% Build the KNN graph to be used later
knn_map_id = zeros(num_data, K); %% KNN label graph
knn_map_dis = zeros(num_data, K); %% KNN distance graph

for i = 1 : num_data
    dataTmp = data(i, :);
    dis = sqrt( sum((ones(num_data, 1) * dataTmp  - data).^2, 2) );
    [c,~] = find(dis == 0);
    dis(c, 1) = 999999;
    [sortv, sortix] = sort(dis, 'ascend');
    knn_map_id(i, :) = sortix(1:K);
    knn_map_dis(i, :) = sortv(1:K);
end


[rdos_score]=knn_cpp(data,knn_map_id,K);

end