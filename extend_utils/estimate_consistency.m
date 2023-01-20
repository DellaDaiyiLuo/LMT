function consist_seg = estimate_consistency(x_measure,x_knnbase,tgrid,k,len)
% x_knnbase = {result_la.xxsamp(1:400,:),result_la.xxsamp(401:800,:)};
% k=6;
% x_measure = result_ori.xxsamp;
d = [0,find(diff(tgrid')>1),numel(tgrid)];
n_seg = numel(d)-1;

n_manifold = numel(x_knnbase);
consist_seg = zeros(n_manifold,n_seg);

for m = 1:n_manifold
    x_tc = x_knnbase{m};
%     display(['k=' num2str(k) ' when searching knn in manifold' num2str(m)]);
    [Idx_m,D] = knnsearch(x_tc,x_measure,'K',k,'Distance','euclidean');
    for j=1:n_seg
        tmpknn_pool = zeros(1,size(x_tc,1));
        for i = d(j)+1:d(j+1)
            tmpknn_pool(Idx_m(i,:)) = max(tmpknn_pool(Idx_m(i,:)),len./(len+D(i,:)));
        end
        consist_seg(m,j) = sum(tmpknn_pool)/(d(j+1)-d(j));
    end
end
end