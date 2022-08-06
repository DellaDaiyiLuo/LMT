function [idx_cong1,idx_cong2,congruence_ori,segcong_ori,m_sc_sm,ci_sc_sm,color] = get_direction_congruence(x_measure,x_knnbase,k,speed_measure,speed_tc,d)
% speed_measure = speed_measure./vecnorm(speed_measure,2,2);
% speed_tc = speed_tc./vecnorm(speed_tc,2,2);

% scale = std(x_knnbase,[],1);
% x_knnbase = x_knnbase./scale;
% x_measure = x_measure./scale;

[Idx,~] = knnsearch(x_knnbase,x_measure,'K',k,'Distance','euclidean');
[Idx_shuf_pool,~] = knnsearch(x_knnbase,x_knnbase,'K',k,'Distance','euclidean');
% figure;hold on
% scatter(x_knnbase(:,1),x_knnbase(:,2),3,'MarkerEdgeColor',[.8,.8,.8])
% scatter(x_measure(:,1),x_measure(:,2),3,'MarkerEdgeColor',[.8,.4,.4])
% scatter(x_measure(165,1),x_measure(165,2),'*')
% scatter(x_knnbase(Idx(165,:),1),x_knnbase(Idx(165,:),2),'o')

knnspeed_mean = zeros(size(speed_measure));
for i=1:size(speed_measure,1)
    knnspeed_mean(i,:)=mean(speed_tc(Idx(i,:),:),1);
end

knnspeed_mean_shuf_pool = zeros(size(speed_tc));
for i=1:size(speed_tc,1)
    knnspeed_mean_shuf_pool(i,:)=mean(speed_tc(Idx_shuf_pool(i,:),:),1);
end

congruence_ori=diag(speed_measure*knnspeed_mean');

segcong_ori=[];
for i=1:numel(d)-1
    segcong_ori = [segcong_ori;mean(congruence_ori(d(i)+1:d(i+1)))];
end

%% shuffled
length = diff(d);
if numel(length)<2
    length_lst=length;
else
    length_lst = [max(min(length)-1,1):max(length)+1]';
end
l = numel(length_lst);
n_sample = sum(length_lst);
n = size(speed_measure,1);
n_tc = size(speed_tc,1);

d_scshuf = [0;cumsum(length_lst)];
segcong_shf_lst = [];
niter=1000;
for iters=1:niter
    ridx1 = randsample(1:n,n_sample,true);
    ridx2 = randsample(1:n_tc,n_sample,true);
    congruence=diag(speed_measure(ridx1,:)*knnspeed_mean_shuf_pool(ridx2,:)');
    segcong=zeros(l,1);
    for i=1:l
        segcong(i) = mean(congruence(d_scshuf(i)+1:d_scshuf(i+1)));
    end
    segcong_shf_lst = [segcong_shf_lst,segcong];
end
ci_sc = zeros(l,2);
m_sc = zeros(l,1);
for i=1:l
    B=sort(segcong_shf_lst(i,:));
    ci_sc(i,1) = B(round(numel(B)*0.975));
    ci_sc(i,2) = B(round(numel(B)*0.025));
    m_sc(i) = mean(B);
end
if numel(length)<2
    idx_cong1=[];idx_cong2=[];
    figure;hold on;
    h=histogram(segcong_shf_lst,30,'Normalization','probability','FaceColor',[0.5,0.5,.5],'EdgeColor','None');
    if segcong_ori>ci_sc(1)
        color=[.7,.2,.2];
    elseif segcong_ori<ci_sc(2)
        color=[.2,.5,.7];
    else
        color=[.3,.3,.3];
    end
    p1=plot([segcong_ori,segcong_ori],ylim,'Color',color);
    plot([ci_sc(1),ci_sc(1)],ylim,'--','Color',[.3,.3,.3])
    plot([ci_sc(2),ci_sc(2)],ylim,'--','Color',[.3,.3,.3])
    xlabel('direction congruence')
    ylabel('relative frequency')
    xl=xlim;
    annotation('doublearrow',(ci_sc-xl(1))/diff(xl),[0.5,0.5])
    annotation('textbox','String','95%','EdgeColor','none','Position',[0.5,0.5,0.1,0.1])
    legend([h,p1],{'random trajectories','PBE example'})
    m_sc_sm = m_sc;
    ci_sc_sm = ci_sc;
else
    color = [];
    m_sc_sm = smoothdata(m_sc,'gaussian',10);
    ci_sc_sm = smoothdata(ci_sc,1,'gaussian',10);
    % lll = repmat(length_lst,1,size(segcong_shf_lst,2));
    figure;hold on
    patch([length_lst;length_lst(end:-1:1)],[ci_sc_sm(:,1);ci_sc_sm(end:-1:1,2)],[.7,.7,.7],'EdgeColor','none','FaceAlpha',.3) 
    plot(length_lst,m_sc_sm,'--','Color',[.5,.5,.5]);
    % scatter(lll(:),segcong_shf_lst(:),'.','MarkerEdgeColor',[.6,.6,.6]);
    idx_cong1 = segcong_ori>ci_sc_sm(length-length_lst(1)+1,1);
    scatter(length(idx_cong1),segcong_ori(idx_cong1),'+','MarkerEdgeColor',[.7,.2,.2]) % significant congruent with tc
    idx_cong2 = segcong_ori<ci_sc_sm(length-length_lst(1)+1,2);
    scatter(length(idx_cong2),segcong_ori(idx_cong2),'+','MarkerEdgeColor',[.2,.5,.7]) % significant in reverse direction with tc
    idx_cong = (~idx_cong1)&(~idx_cong2);
    scatter(length(idx_cong),segcong_ori(idx_cong),'+','MarkerEdgeColor',[.3,.3,.3]) % not significant
    xlim([length_lst(1),length_lst(end)])
    % plot(length_lst,ci_sc_sm(:,1)');
    % plot(length_lst,ci_sc_sm(:,2)');
    ylabel('direction congruence')
    xlabel('number of steps of trajectories')
    legend({'95% interval','mean of shuffled','same direction','reverse direction','insignificant'})


    figure;
    hold on
    c = diff(x_measure(:,1:2));
    scatter(x_knnbase(:,1),x_knnbase(:,2),3,'MarkerEdgeColor',[.8,.8,.8]);
    for i=1:numel(d)-1
        if idx_cong1(i)
            p1=quiver(x_measure(d(i)+1:d(i+1)-1,1),x_measure(d(i)+1:d(i+1)-1,2),c(d(i)+1:d(i+1)-1,1),c(d(i)+1:d(i+1)-1,2),0,'LineWidth',0.5,'Color',[.7,.2,.2]);
        elseif idx_cong2(i)
            p2=quiver(x_measure(d(i)+1:d(i+1)-1,1),x_measure(d(i)+1:d(i+1)-1,2),c(d(i)+1:d(i+1)-1,1),c(d(i)+1:d(i+1)-1,2),0,'LineWidth',0.5,'Color',[.2,.5,.7]);
        else
            p3=quiver(x_measure(d(i)+1:d(i+1)-1,1),x_measure(d(i)+1:d(i+1)-1,2),c(d(i)+1:d(i+1)-1,1),c(d(i)+1:d(i+1)-1,2),0,'LineWidth',0.5,'Color',[.3,.3,.3]);
        end
    end

%     legend([p1,p2,p3],{'same direction','reverse direction','insignificant'})

    % %% shuffled
    % segcong_lst = [];
    % n = numel(d)-1;
    % l=size(speed_measure,1);
    % niter=1000;
    % for iters=1:niter
    %     ridx = randperm(l);
    %     congruence=diag(speed_measure(ridx,:)*knnspeed_mean');
    %     segcong=zeros(n,1);
    %     for i=1:n
    %         segcong(i) = mean(congruence(d(i)+1:d(i+1)));
    %     end
    %     segcong_lst = [segcong_lst,segcong];
    % end
    % 
    % 
    % % plot
    % segcong_pool=segcong_lst(:);
    % B = sort(segcong_pool);
    % bar1 = B(round(numel(B)*0.975));
    % bar2 = B(round(numel(B)*0.025));
    % figure;h=histogram(segcong_pool,50,'Normalization','probability','FaceColor',[0.5,0.5,.5],'EdgeColor','None');
    % yl=ylim;
    % hold on
    % for i=1:numel(segcong_ori)
    %     if segcong_ori(i)>0
    %         p1=plot([segcong_ori(i) segcong_ori(i)],yl,'Color',[.7,.2,.2],'LineWidth',1);
    %     else
    %         p2=plot([segcong_ori(i) segcong_ori(i)],yl,'Color',[.2,.5,.7],'LineWidth',1);
    %     end
    % end
    % p3 = plot([bar1 bar1],yl,'k--','LineWidth',2);
    % p3 = plot([bar2 bar2],yl,'k--','LineWidth',2);
    % ratio = [numel(find(segcong_ori>0)),numel(segcong_ori)];
    % if diff(ratio)==0
    %     legend(p1,{'segment congruence>0'})
    % else
    %     legend([h,p1,p2],{'shuffled input data','original segcong>0','original segcong<=0'})
    % end
    % title('Histogram of segment congruences')


    % %%
    % figure;
    % hold on
    % if numel(segcong_ori)>1
    %     color=(segcong_ori-min(segcong_ori))/(max(segcong_ori)-min(segcong_ori));
    % else
    %     color=1;
    % end
    % c = diff(x_measure(:,1:2));
    % scatter(x_knnbase(:,1),x_knnbase(:,2),3,'MarkerEdgeColor',[.8,.8,.8]);
    % for i=1:numel(d)-1
    %     if segcong_ori(i)>1.5
    %         quiver(x_measure(d(i)+1:d(i+1)-1,1),x_measure(d(i)+1:d(i+1)-1,2),c(d(i)+1:d(i+1)-1,1),c(d(i)+1:d(i+1)-1,2),0,'LineWidth',0.5,'Color',[1,0.5,0.5]*color(i))
    %     end
    % end
    % colorbar('Ticks',[0,1],...
    %          'TickLabels',{num2str(min(segcong_ori)),num2str(max(segcong_ori))})
end
end
