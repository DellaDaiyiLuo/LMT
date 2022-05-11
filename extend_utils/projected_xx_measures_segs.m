function [stepdis_seg, knndis_seg, knnportion_seg] = projected_xx_measures_segs(x_measure,x_knnbase,d,plot_result,k)
% Measures for unsupervised learned result
% x_measure: projected xx to be measured
% x_knnbase: tuning curve xx to be the base for x_measure searching knn
% d: point indices when time is not continous, corresponding path segments
% will be removed from path length
% plot_result: plot result or not

disdiff = vecnorm(diff(x_measure),2,2);
n_seg = numel(d)-1;
stepdis_seg = zeros(n_seg,1);
for i=1:n_seg
    segidx = d(i)+1:d(i+1)-1;
    stepdis_seg(i) = sum(disdiff(segidx))/numel(segidx);
end
n_manifold = numel(x_knnbase);

knndis_seg = zeros(n_manifold,n_seg);
knnportion_seg = zeros(n_manifold,n_seg);
knn_in_tc = {};
for m = 1:n_manifold
    x_tc = x_knnbase{m};
    display(['k=' num2str(k) ' when searching knn in manifold' num2str(m)]);
    [Idx,D] = knnsearch(x_tc,x_measure,'K',k,'Distance','euclidean');
    for i=1:n_seg
        segidx = d(i)+1:d(i+1)-1;
        knndis_seg(m,i) = mean(D(segidx,:),'all'); % knn distances to tc manifold
        tmpidx = Idx(segidx,:);
        tmptbl = tabulate(tmpidx(:));
        knnportion_seg(m,i) = sum(tmptbl(:,2)>0)/size(x_tc,1); % knn distribution on tc
    end
    knn_in_tc = [knn_in_tc, {tmptbl(:,2)>0}];
end

if plot_result
    figure;hold on
    p = [];
    labels = {};
    for m=1:n_manifold
        p1=scatter(x_knnbase{m}(:,1),x_knnbase{m}(:,2),7,'MarkerFaceColor',[.8,.8,.8]*m/n_manifold,'MarkerEdgeColor','None');
        p = [p,p1]; 
        labels = [labels,{['TC manifold', num2str(m)]}];
        p2=scatter(x_knnbase{m}(knn_in_tc{m},1),x_knnbase{m}(knn_in_tc{m},2),10,'*');
        p = [p,p2]; 
        labels = [labels,{['knn on manifold', num2str(m)]}];
    end
    p3=plot(x_measure(d(i)+1:d(i+1),1),x_measure(d(i)+1:d(i+1),2),'Color',[0.80,0.24,0.80]);%[0.7,0.5,0.5]); %
    
    legend([p,p3],[labels,{'projected data'}])
end
end