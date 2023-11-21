function measures_hist_multiple(id,cellknndis_lst,timeknndis_lst,knndis_ori,name)
hold on;
limit = [floor(min([cellknndis_lst{1}(:,id);timeknndis_lst{1}(:,id)])),ceil(max([cellknndis_lst{1}(:,id);timeknndis_lst{1}(:,id)]))];
histogram(cellknndis_lst{1}(:,id),25,'Normalization','probability','BinLimits',limit,'DisplayStyle','stairs','EdgeColor',[.1,.4,.6])
histogram(timeknndis_lst{1}(:,id),25,'Normalization','probability','BinLimits',limit,'FaceColor',[.1,.4,.6],'EdgeColor','None')
yl = ylim;
plot([knndis_ori(1,id) knndis_ori(1,id)],yl,'--','Color',[.1,.4,.6]*.5,'LineWidth',2);

limit = [floor(min([cellknndis_lst{2}(:,id);timeknndis_lst{2}(:,id)])),ceil(max([cellknndis_lst{2}(:,id);timeknndis_lst{2}(:,id)]))];
histogram(cellknndis_lst{2}(:,id),25,'Normalization','probability','BinLimits',limit,'DisplayStyle','stairs','EdgeColor',[.7,0.1,0.2])
histogram(timeknndis_lst{2}(:,id),25,'Normalization','probability','BinLimits',limit,'FaceColor',[.7,0.1,0.2],'EdgeColor','None')
plot([knndis_ori(2,id) knndis_ori(2,id)],yl,'--','Color',[.7,0.1,0.2]*.5,'LineWidth',2);
% scatter(knndis_ori(1),yl(2)*0.97,100,'p','filled','MarkerFaceColor', [0.5,0.15,.5])
legend({'M1 cell-shuf','      time-shuf','      original','M2 cell-shuf','      time-shuf','      original','original context'})
xlabel(name)
ylabel('relative frequency')
end