function [speed,d] = get_speed(tgrid,xplot,plot_quiver)
d = [0 find(diff(tgrid')>1) numel(tgrid)];
speed = [];
% for i=1:numel(d)-1
%         x_tmp = xplot(d(i)+1:d(i+1),:);
%         x_sm = smoothdata(x_tmp,1,'gaussian',10);%20230103 %10);%x_tmp;%
%         c = diff(x_sm);
%         speed = [speed;c;c(end,:)];
% end
for i=1:numel(d)-1
        x_tmp = xplot(d(i)+1:d(i+1),:);
        s = gradient(x_tmp');
%         x_sm = smoothdata(x_tmp,'gaussian',10);
%         s = gradient(x_sm');
        speed = [speed;s'];
end
% for dim=1:size(xplot,2)
%     s=[];
%     for i=1:numel(d)-1
%         x_tmp = xplot(d(i)+1:d(i+1),dim);
%         c = gradient(x_tmp);
% %         x_sm = smoothdata(x_tmp,'gaussian',5);
% %         c = gradient(x_sm);
%         s = [s;c];
%     end
%     speed = [speed,s];
% end

if plot_quiver
    switch size(xplot,2)
        case 2
            figure;quiver(xplot(:,1),xplot(:,2),speed(:,1),speed(:,2))
        case 3
            quiver3(xplot(:,1),xplot(:,2),xplot(:,3),speed(:,1),speed(:,2),speed(:,3),3)
    end
end
end