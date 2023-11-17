function [dircong_seg,spdcong_seg] = get_dir_congruence(tgrid,x_measure,x_knnbase_all,speed_tc,k)

[speed_measure,d] = get_speed(tgrid,x_measure,0);
disdiff = vecnorm(speed_measure,2,2);

[Idx,D] = knnsearch(x_knnbase_all,x_measure,'K',k,'Distance','euclidean');
w = 1./(1+D);

congruence = zeros(size(speed_measure,1),1); % only consider line
congruence_spd = zeros(size(speed_measure,1),1); % consider same/opposite direction
for i=1:size(speed_measure,1)
    congruence(i)=w(i,:)*abs(speed_tc(Idx(i,:),:)*speed_measure(i,:)')/disdiff(i);
    congruence_spd(i)=w(i,:)*(speed_tc(Idx(i,:),:)*speed_measure(i,:)')/disdiff(i);
end

dircong_seg=[];
spdcong_seg=[];
for i=1:numel(d)-1
    segidx = d(i)+1:d(i+1);
    dircong_seg = [dircong_seg,mean(congruence(segidx))];
    spdcong_seg = [spdcong_seg,mean(congruence_spd(segidx))];
end
end