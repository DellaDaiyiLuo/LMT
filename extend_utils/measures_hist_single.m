function measures_hist_single(id,cell_lst,time_lst,ori,name)
figure;
hold on;
limit = [floor(min([cell_lst(:,id);time_lst(:,id)])),ceil(max([cell_lst(:,id);time_lst(:,id)]))];
histogram(cell_lst(:,id),25,'Normalization','probability','BinLimits',limit,'DisplayStyle','stairs','EdgeColor',[0.5,0.15,.5])
histogram(time_lst(:,id),25,'Normalization','probability','BinLimits',limit,'FaceColor',[0.5,0.15,.5],'EdgeColor','None')
plot([ori(id) ori(id)],ylim,'--','Color',[0.5,0.15,.5]*0.7,'LineWidth',2);
xlabel(name)
ylabel('relative frequency')
legend({'cell-shuf','time-shuf','original'})
end