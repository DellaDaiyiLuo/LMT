addpath(genpath('/home/kemerelab/Downloads/EEG/LMT-master'));warning off
addpath(genpath('/home/kemerelab/Downloads/EEG/extend_utils'));

% -------- load tc data ----------- %
% run cross valid
load('/home/kemerelab/Downloads/EEG/Elife2018_data/result_elife_t1_1000_2000_nf2_init_pca.mat')
xbackground={result_la.xxsamp};
% run mix tracks
load('/home/kemerelab/Downloads/EEG/Elife2018_data/True_Shuffle_MIX/result_elife_mix_nf2_init_pca.mat')
load('/home/kemerelab/Downloads/EEG/Elife2018_data/True_Shuffle_MIX/result_elife_mix_nf2_init_pca_xbackground.mat')


tgrid_tc = setopt.tgrid;
xx_tc = xx;

% -------- load test data ---------%
% simul
spikes = spikes(:,1:60)';

% track1 run
load('/home/kemerelab/Downloads/EEG/Elife2018_data/track1_100ms.mat')
position=position(idx);
spikes=spikes(:,idx);
% PBE
load('/home/kemerelab/Downloads/EEG/Elife2018_data/track1_PBEs_4ms.mat')
load('/home/kemerelab/Downloads/EEG/figures/PBE/pberatio.mat')
position=xinit_ratio';
% -------- define input data ----------- %


portion = [2001:2500];%[4739:4793];%PBE2  %%[4296:4350];%PBE1;

xx0 = position(portion)'; 
yy0=double(spikes(:,portion))';
tgrid = double(idx(portion)');%portion'; %
unique(diff(tgrid)) % pbe segment
%% ---------- run projection -----------%%
niter = 20;
shuffle = 'no'; 
[result, setopt, fftc, stepdis_ori, knndis_ori, knnportion_ori,xx,order] = run_projection(xx0,yy0,...
    tgrid,result_la,niter,'shuffle',shuffle,'xbackground',xbackground,'k',k,'tctype','run','scaler',scaler,'newdatatype','pbe','segmeasure',0); % two tracks pbe

% [result, setopt, fftc, stepdis_ori, knndis_ori, knnportion_ori,xx,order] = run_projection(xx0,yy0,tgrid,result_la,niter,'shuffle',shuffle,'xbackground',xbackground,'k',k); % two tracks
% [result, setopt, fftc, pathlen, knndis, knnportion,xx,order] =
% run_projection(xx0,yy0,tgrid,result_la,niter,'shuffle',shuffle); % run
% cross validation

%----find PBE we want ------%


figure;
d = [0 find(diff(tgrid')>1) numel(tgrid)];
lllll = diff(d)-1;
panel=[knnportion_ori(1,:)./lllll*1000;lllll]; % knnportion/number of steps

i=1;

seg=d(i)+1:d(i+1);
clf
scatter(result_la.xxsamp(:,1),result_la.xxsamp(:,2),3,[.5,.5,.5])
hold on
plot(result.xxsamp(seg,1),result.xxsamp(seg,2));
scatter(setopt.xplds(seg,1),setopt.xplds(seg,2),3,xinit_ratio(seg));
title(i)
numel(seg)
portion([seg(1),seg(end)])
i=i+1;
%---------------------------%

figure;
show_latent_variable(setopt.xplds,xx,xbackground,tgrid,'scatter_only',0)
title({'Bayesian initialization of new input data'})
axis([-50 50 -55 55])
axis([-83,63,-55,50])

figure;
show_latent_variable(result.xxsamp,xx,xbackground,tgrid)
title({'New input data projected into latent space'})
axis([-50 50 -55 55])
axis([-83,63,-55,50])

save('decodedPBE2.mat', 'result', 'portion', 'setopt','xx','order','stepdis_ori', 'knndis_ori', 'knnportion_ori');

% -------- additional projection runs------- %
niter = 10;
xxtc = result_la.xxsamp;
[result_20,~] = additional_projection(xx0,yy0,fftc,xxtc,setopt,result,niter);
figure;
show_latent_variable(result_20.xxsamp,xx,xbackground,tgrid)
title('projected x, iter = 20')



%% ----------- subsequent analysis ------------%
% ---------- direction -----------%
tgrid = setopt.tgrid;
x_measure = result.xxsamp;
[speed_measure,d] = get_speed(tgrid,x_measure,1);
x_knnbase = result_la.xxsamp;
[speed_tc, d_tc] = get_speed(tgrid_tc,x_knnbase,1);
k=6;
% for PBE, put in single PBE is better
[congruence_ori,segcong_ori,~,~,color] = get_direction_congruence(x_measure,x_knnbase,k,speed_measure,speed_tc,d);

% [congruence_ori,segcong_ori,m_sc_sm,ci_sc_sm,color] = get_direction_congruence(x_measure,x_knnbase,k,speed_measure,speed_tc,d);


% ------------ PBE plotting ------------------ %
open('/home/kemerelab/Downloads/EEG/figures/TTF/TTF2.fig');
hold on
% show_latent_variable(result.xxsamp(1:32,:),xx(1:32),[],tgrid,'line_only',1,'line_color',[.7,.2,.2])
% show_latent_variable(result.xxsamp(33:end,:),xx(33:end),[],tgrid,'line_only',1,'line_color',[.2,.5,.7])
% show_latent_variable(setopt.xplds,xx,[],tgrid,'line_only',1,'line_color',[.7,.7,.7])
show_latent_variable(result.xxsamp,xx,[],setopt.tgrid,'line_only',1)
interval = 4;
% quiver(result.xxsamp(1:interval:end,1),result.xxsamp(1:interval:end,2),speed_measure(1:interval:end,1),speed_measure(1:interval:end,2),'Color',color)
quiver(result.xxsamp(1:interval:32,1),result.xxsamp(1:interval:32,2),speed_measure(1:interval:32,1),speed_measure(1:interval:32,2),'Color',[.7,.2,.2])
quiver(result.xxsamp(33:interval:end,1),result.xxsamp(33:interval:end,2),speed_measure(33:interval:end,1),speed_measure(33:interval:end,2),'Color',[.2,.5,.7])
title('PBE example')

figure
scatter(1:numel(tgrid),xinit_ratio(portion),15,xinit_ratio(portion),'filled')
axis([0,numel(tgrid)+1,-1,266])
ylabel('Bayesian decoded position')
xlabel('time bin')
title('PBE example')

% figure;hold on;histogram(congruence,'FaceColor',[.5,.5,.5],'EdgeColor','None');
% xlabel('speed congruence')
% ylabel('count')
% title('Histogram of speed congruence')
% l = 400;
% for i=1:numel(segcong)
%     if segcong(i)>0
%         p1=plot([segcong(i) segcong(i)],[0 l],'Color',[.7,.2,.2],'LineWidth',1);
%     else
%         p2=plot([segcong(i) segcong(i)],[0 l],'Color',[.2,.5,.7],'LineWidth',1);
%     end
% end
% p3 = plot([0 0],[0 l],'k--','LineWidth',2);
% ratio = [numel(find(segcong>0)),numel(segcong)]
% if diff(ratio)==0
%     legend([p1,p3],{'segment congruence>0','threshold'})
% else
%     legend([p1,p2,p3],{'segment congruence>0','segment congruence<=0','threshold'})
% end

% ---------- shuffled data pool ------------- %
n_manifold = numel(xbackground);
stepdis_lst = [];
knndis_lst = [];
knnportion_lst = [];
% for i=1:n_manifold
%     knndis_lst{i}=[];
%     knnportion_lst{i}=[];
% end

niter = 10;
shuffle = 'cell'; 
load([shuffle 'shuffle_pool.mat'])
pool_size = 100;
tic
for ff=1:5
    for times = 1:pool_size
        % pbe
        [~, ~, ~, stepdis, knndis, knnportion,~,~] = run_projection(xx0,yy0,...
            tgrid,result_la,niter,'shuffle',shuffle,'draw',0,'xbackground',xbackground,'k',k,'tctype','run','scaler',scaler,'newdatatype','pbe'); % two tracks pbe
%         stepdis_lst = [stepdis_lst, stepdis];
%         for i=1:n_manifold
%             knndis_lst{i} = [knndis_lst{i}; knndis(i,:)];
%             knnportion_lst{i} = [knnportion_lst{i}; knnportion(i,:)];
%         end
        % run
%         [~, ~, ~, stepdis, knndis, knnportion,~,~] = run_projection(xx0,yy0,tgrid,result_la,niter,'shuffle',shuffle,'draw',1,'xbackground',xbackground,'k',k);
        stepdis_lst = [stepdis_lst, stepdis];
        knndis_lst = [knndis_lst, knndis];
        knnportion_lst = [knnportion_lst, knnportion];
    end
    save([shuffle, 'shuffle_pool_pbe2.mat'],'stepdis_lst', 'knndis_lst', 'knnportion_lst');
end
toc

timestepdis_lst=stepdis_lst;
timeknndis_lst=knndis_lst;
timeknnportion_lst=knnportion_lst;
load('cellshuffle_pool.mat')

measures_hist(stepdis_lst,timestepdis_lst,stepdis_ori,knndis_lst,timeknndis_lst,knndis_ori,knnportion_lst,timeknnportion_lst,knnportion_ori)
title('PBE example1')

timeknndisratio=timeknndis_lst(1,:)./timeknndis_lst(2,:);
cellknndisratio=knndis_lst(1,:)./knndis_lst(2,:);
oriknndisratio=knndis_ori(1,:)./knndis_ori(2,:);
figure;
hold on;
limit = [floor(min([cellknndisratio,timeknndisratio])),ceil(max([cellknndisratio,timeknndisratio]))];
histogram(cellknndisratio,100,'Normalization','probability','BinLimits',limit,'DisplayStyle','stairs','EdgeColor',[0.5,0.15,.5])
histogram(timeknndisratio,100,'Normalization','probability','BinLimits',limit,'FaceColor',[0.5,0.15,.5],'EdgeColor','None')
plot([oriknndisratio oriknndisratio],ylim,'--','Color',[0.5,0.15,.5]*0.7,'LineWidth',2);
xlabel('average step distance')
ylabel('relative frequency')
legend({'cell-shuf','time-shuf','original'})






% -------------- evolution -----------------%
[pks,locs] = findpeaks(xx,'MinPeakHeight',0.6*max(xx),'MinPeakDistance', 80);
segs = [1;locs;numel(xx)];
n = numel(segs)-1;
figure;plot(xx,'k');hold on;
for i=1:n-1
    plot([locs(i),locs(i)],[0 270],'r--')
end
axis([0,2785,0,270])

[speed_lin,~] = get_speed(tgrid,xx,0);
m = max(abs(speed_lin));
xplot = result_la.xxsamp; %result.xxsamp;
figure
for i=1:n-1
    subplot(4,5,i)
    title(['bout' num2str(i)])
    show_latent_variable(xplot,abs(speed_lin),{result_la.xxsamp},tgrid,'part',segs(i):segs(i+1),'legend',0);
    caxis([0,m])
end
i=n;
subplot(4,5,i)
title(['bout' num2str(i)])
show_latent_variable(xplot,abs(speed_lin),{result_la.xxsamp},tgrid,'part',segs(i):segs(i+1),'legend',1);
caxis([0,m])









