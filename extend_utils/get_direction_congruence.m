function [congruence,segcong] = get_direction_congruence(x_measure,x_knnbase,k,speed_measure,speed_tc,d)
% speed_measure = speed_measure./vecnorm(speed_measure,2,2);
% speed_tc = speed_tc./vecnorm(speed_tc,2,2);

% scale = std(x_knnbase,[],1);
% x_knnbase = x_knnbase./scale;
% x_measure = x_measure./scale;

[Idx,~] = knnsearch(x_knnbase,x_measure,'K',k,'Distance','euclidean');

% figure;hold on
% scatter(x_knnbase(:,1),x_knnbase(:,2),3,'MarkerEdgeColor',[.8,.8,.8])
% scatter(x_measure(:,1),x_measure(:,2),3,'MarkerEdgeColor',[.8,.4,.4])
% scatter(x_measure(165,1),x_measure(165,2),'*')
% scatter(x_knnbase(Idx(165,:),1),x_knnbase(Idx(165,:),2),'o')

knnspeed_mean = zeros(size(speed_measure));
for i=1:size(speed_measure,1)
    knnspeed_mean(i,:)=mean(speed_tc(Idx(i,:),:),1);
end
congruence=diag(speed_measure*knnspeed_mean');
figure;plot(congruence)

segcong=[];
for i=1:numel(d)-1
    segcong = [segcong;mean(congruence(d(i)+1:d(i+1)))];
end

figure;
hold on
if numel(segcong)>1
    color=(segcong-min(segcong))/(max(segcong)-min(segcong));
else
    color=1;
end
c = diff(x_measure(:,1:2));
scatter(x_knnbase(:,1),x_knnbase(:,2),3,'MarkerEdgeColor',[.8,.8,.8]);
for i=1:numel(d)-1
    quiver(x_measure(d(i)+1:d(i+1)-1,1),x_measure(d(i)+1:d(i+1)-1,2),c(d(i)+1:d(i+1)-1,1),c(d(i)+1:d(i+1)-1,2),0,'LineWidth',0.5,'Color',[1,0.5,0.5]*color(i))
end
colorbar('Ticks',[0,1],...
         'TickLabels',{num2str(min(segcong)),num2str(max(segcong))})
end
