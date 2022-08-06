function [dircong_seg] = get_dir_congruence(tgrid,x_measure,x_knnbase_all,speed_tc,k)
[speed_measure,d] = get_speed(tgrid,x_measure,0);

[Idx,D] = knnsearch(x_knnbase_all,x_measure,'K',k,'Distance','euclidean');
w = 1./(2+D);

congruence = zeros(size(speed_measure,1),1);
for i=1:size(speed_measure,1)
    congruence(i)=w(i,:)*abs(speed_tc(Idx(i,:),:)*speed_measure(i,:)')/norm(speed_measure(i,:));
end

dircong_seg=[];
for i=1:numel(d)-1
    dircong_seg = [dircong_seg,mean(congruence(d(i)+1:d(i+1)))];
end
end