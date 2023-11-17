function [knndis_seg, uniqueknn_seg,knnportion_seg_m] = projected_xx_measures_segs(x_measure,x_knnbase,d,plot_result,k)
% Measures for unsupervised learned result
% x_measure: projected xx to be measured
% x_knnbase: tuning curve xx to be the base for x_measure searching knn
% d: point indices when time is not continous, corresponding path segments
% will be removed from path length
% plot_result: plot result or not

% % stepdis
% disdiff = vecnorm(diff(x_measure)',2);
n_seg = numel(d)-1;
% stepdis_seg = zeros(1,n_seg);
% for i=1:n_seg
%     segidx = d(i)+1:d(i+1)-1;
%     stepdis_seg(i) = sum(disdiff(segidx))/numel(segidx);
% end

% knndis, knnportion
n_manifold = numel(x_knnbase);
knndis_seg = zeros(n_manifold,n_seg);
knnportion_seg_m = zeros(n_manifold,n_seg);
uniqueknn_seg = zeros(1,n_seg);
knn_in_tc = {};
x_knnbase_all = [];
for i_m=1:n_manifold
    x_knnbase_all = [ x_knnbase_all; x_knnbase{i_m}];
end
[Idx,~] = knnsearch(x_knnbase_all,x_measure,'K',k,'Distance','euclidean');
for m = 1:n_manifold
    x_tc = x_knnbase{m};
    display(['k=' num2str(k) ' when searching knn in manifold' num2str(m)]);
    [Idx_m,D] = knnsearch(x_tc,x_measure,'K',k,'Distance','euclidean');
    for i=1:n_seg
        segidx = d(i)+1:d(i+1);
        knndis_seg(m,i) = mean(D(segidx,:),'all'); % knn distances to tc manifold
        tmpidx = Idx(segidx,:);
        tmptbl = tabulate(tmpidx(:));
        uniqueknn_seg(i) = sum(tmptbl(:,2)>0)/(d(i+1)-d(i));
        tmpidx = Idx_m(segidx,:);
        tmptbl = tabulate(tmpidx(:));
        knnportion_seg_m(m,i) = sum(tmptbl(:,2)>0)/(d(i+1)-d(i)); % knn distribution on tc
%         knnportion_seg_m(m,i) = sum(tmptbl(:,2)>0)/size(x_tc,1); % knn distribution on tc
    end
    knn_in_tc = [knn_in_tc, {tmptbl(:,2)>0}];
end

% % seg_congruence _new
% [dircong_seg,spdcong_seg] = get_dir_congruence(tgrid,x_measure,x_knnbase_all,speed_tc,k);


% % seg_congruence _old
% [speed_measure,d] = get_speed(tgrid,x_measure,0);
% 
% 
% knnspeed_mean = zeros(size(speed_measure));
% for i=1:size(speed_measure,1)
%     knnspeed_mean(i,:)=mean(speed_tc(Idx(i,:),:),1);
% end
% 
% congruence_ori=diag(speed_measure*knnspeed_mean');
% 
% spdcong_seg=[];
% for i=1:numel(d)-1
%     spdcong_seg = [spdcong_seg,mean(congruence_ori(d(i)+1:d(i+1)))];
% end

% plot
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