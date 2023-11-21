%% plot single environment trajectory
figure;hold on;
a1=area([1000,2000],[0,265; 0,265],'FaceColor',[1,1,.9],'EdgeColor','None');
a2=area([2001,2785],[0,265; 0,265],'FaceColor',[1,.95,.85],'EdgeColor','None');
d = [0 find(diff(idx)>1) numel(idx)];
for i=1:numel(d)-1
    p1=plot(d(i)+1:d(i+1),position(d(i)+1:d(i+1)),'Color',[.7,.5,.5]);
end
scatter(1:2785,position,3,position)
axis([0,2785,0,265])
legend([a1(1),a2(1)],{'Chunk1','Chunk2'})

%% plot two environments trajectory
plotst2=501;
idx = [idx1,idx2(plotst2:end)];
position = [position1,position2(plotst2:end)];
figure;ax1 = axes;
hold on
a1=area([portion1(1),portion1(end)],[0,265; 0,265],'FaceColor',[1,1,.9],'EdgeColor','None');
area([portion2(1),portion2(end)]+numel(idx1)-plotst2,[0,265; 0,265],'FaceColor',[1,1,.9],'EdgeColor','None');
a2=area([2001,2785],[0,265; 0,265],'FaceColor',[.9,.95,1],'EdgeColor','None');
d = [0 find(diff(idx1)>1) numel(idx1)];
for i=1:numel(d)-1
    p1=plot(d(i)+1:d(i+1),position(d(i)+1:d(i+1)),'k');
end
plot([numel(idx1),numel(idx1)],[0,265],'--','Color',[.6,.6,.6])
d = [0 find(diff(idx2(plotst2:end))>1) numel(idx2(plotst2:end))]+numel(idx1);
for i=1:numel(d)-1
    p2=plot(d(i)+1:d(i+1),position(d(i)+1:d(i+1)),'Color',[.6,.6,.6]);
end
% ax1 = axes;
scatter(ax1,portion1,xx1,10,xx1,'filled')
scatter(ax1,2001:2785,position1(2001:2785),10,position1(2001:2785),'filled')
xlabel('time bin')
ylabel('animal position')
ax2 = axes;
scatter(ax2,portion2+numel(idx1)-plotst2,xx2,10,xx2,'filled')
linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'parula')
colormap(ax2,'copper')
axis([1200,4000,0,265])


set([ax1,ax2],'Position',[.19 .2 .69 .715]);
cb1 = colorbar(ax1,'Position',[.1 .11 0.036 .815]);
cb1.Label.String='track1 position';
cb2 = colorbar(ax2,'Position',[.9 .11 0.036 .815]);
cb2.Label.String='track2 position';
legend([a1(1),a2(1)],{'train','test'})








