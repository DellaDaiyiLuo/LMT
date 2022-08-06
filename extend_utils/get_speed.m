function [speed,d] = get_speed(tgrid,xplot,plot_quiver)
d = [0 find(diff(tgrid')>1) numel(tgrid)];
speed = [];
for i=1:numel(d)-1
    x_tmp = xplot(d(i)+1:d(i+1),:);
    x_sm = smoothdata(x_tmp,1,'gaussian',10);
    c = diff(x_sm);
    speed = [speed;c;c(end,:)];
end
if plot_quiver
    switch size(xplot,2)
        case 2
            figure;quiver(xplot(:,1),xplot(:,2),speed(:,1),speed(:,2))
        case 3
            quiver3(xplot(:,1),xplot(:,2),xplot(:,3),speed(:,1),speed(:,2),speed(:,3),3)
    end
end
end