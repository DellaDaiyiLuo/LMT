function measures_hist(cellstepdis_lst,timestepdis_lst,stepdis_ori,cellknndis_lst,timeknndis_lst,knndis_ori,cellknnportion_lst,timeknnportion_lst,knnportion_ori)
figure;
hold on;
limit = [floor(min([cellstepdis_lst,timestepdis_lst])),ceil(max([cellstepdis_lst,timestepdis_lst]))];
histogram(cellstepdis_lst,25,'Normalization','probability','BinLimits',limit,'DisplayStyle','stairs','EdgeColor',[0.5,0.15,.5])
histogram(timestepdis_lst,25,'Normalization','probability','BinLimits',limit,'FaceColor',[0.5,0.15,.5],'EdgeColor','None')
plot([stepdis_ori stepdis_ori],ylim,'--','Color',[0.5,0.15,.5]*0.7,'LineWidth',2);
xlabel('average step distance')
ylabel('relative frequency')
legend({'cell-shuf','time-shuf','original'})

figure;
hold on;
limit = [floor(min([cellknndis_lst(:);timeknndis_lst(:)])),ceil(max([cellknndis_lst(:);timeknndis_lst(:)]))];
histogram(cellknndis_lst(1,:),25,'Normalization','probability','BinLimits',limit,'DisplayStyle','stairs','EdgeColor',[.1,.4,.6])
histogram(timeknndis_lst(1,:),25,'Normalization','probability','BinLimits',limit,'FaceColor',[.1,.4,.6],'EdgeColor','None')
yl = ylim;
plot([knndis_ori(1) knndis_ori(1)],yl,'--','Color',[.1,.4,.6]*.5,'LineWidth',2);

histogram(cellknndis_lst(2,:),25,'Normalization','probability','BinLimits',limit,'DisplayStyle','stairs','EdgeColor',[.7,0.1,0.2])
histogram(timeknndis_lst(2,:),25,'Normalization','probability','BinLimits',limit,'FaceColor',[.7,0.1,0.2],'EdgeColor','None')
plot([knndis_ori(2) knndis_ori(2)],yl,'--','Color',[.7,0.1,0.2]*.5,'LineWidth',2);
scatter(knndis_ori(1),yl(2)*0.97,100,'p','filled','MarkerFaceColor', [0.5,0.15,.5])
legend({'M1 cell-shuf','      time-shuf','      original','M2 cell-shuf','      time-shuf','      original','original context'})
xlabel('average knn distance')
ylabel('relative frequency')

figure;
hold on;
limit = [min([cellknnportion_lst(:);timeknnportion_lst(:)]),max([cellknnportion_lst(:);timeknnportion_lst(:)])];
histogram(cellknnportion_lst(1,:),25,'Normalization','probability','BinLimits',limit,'DisplayStyle','stairs','EdgeColor',[.1,.4,.6])
histogram(timeknnportion_lst(1,:),25,'Normalization','probability','BinLimits',limit,'FaceColor',[.1,.4,.6],'EdgeColor','None')
yl = [0,.9];
plot([knnportion_ori(1) knnportion_ori(1)],yl,'--','Color',[.1,.4,.6]*.5,'LineWidth',2);

histogram(cellknnportion_lst(2,:),25,'Normalization','probability','BinLimits',limit,'DisplayStyle','stairs','EdgeColor',[.7,0.1,0.2])
histogram(timeknnportion_lst(2,:),25,'Normalization','probability','BinLimits',limit,'FaceColor',[.7,0.1,0.2],'EdgeColor','None')
plot([knnportion_ori(2) knnportion_ori(2)],yl,'--','Color',[.7,0.1,0.2]*.5,'LineWidth',2);
scatter(knnportion_ori(1),yl(2)*0.97,100,'p','filled','MarkerFaceColor', [0.5,0.15,.5])
% legend({'track1 cell-shuf','           time-shuf','           original','track2 cell-shuf','           time-shuf','           original'})

xlabel('knn portion')
ylabel('relative frequency')
end
