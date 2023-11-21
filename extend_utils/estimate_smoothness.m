function [smoothness,stepdis_seg] = estimate_smoothness(xplot,d,norm_dis)
n_seg = numel(d)-1;
smoothness = zeros(1,n_seg);
if norm_dis
    for i=1:n_seg
        x = xplot(d(i)+1:d(i+1),:);
        dx = diff(x);
        T = dx./vecnorm(dx,2,2);
        dT = diff(T);
        smoothness(i) = -mean(vecnorm(dT,2,2)./vecnorm(dx(1:end-1,:),2,2)); % negative curvature
    end
else
    for i=1:n_seg
        x = xplot(d(i)+1:d(i+1),:);
        smoothness(i) = -mean(vecnorm(diff(x,2),2,2));
    end
end

disdiff = vecnorm(diff(xplot)',2);
stepdis_seg = zeros(1,n_seg);
for i=1:n_seg
    segidx = d(i)+1:d(i+1)-1;
    stepdis_seg(i) = sum(disdiff(segidx))/numel(segidx);
end
end