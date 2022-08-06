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
tb = 12;
load(['/home/kemerelab/Downloads/EEG/Elife2018_data/track1_PBEs_' num2str(tb) 'ms.mat'])
scaler=tb/100;

load('/home/kemerelab/Downloads/EEG/figures/PBE/pberatio.mat')
position=xinit_ratio';
portion = [4739:4793];%PBE2  %%[4296:4350];%PBE1;[event_edge_sw(501,1):event_edge_sw(560,2)];  %[10000:17526];%

% -------- define input data ----------- %
% for maze-----
maze = 1;
load('result_la_maze12_500ms_3nf_30iter_allrun_idxtgrid.mat')
xbackground={result_la.xxsamp(1:815,:),result_la.xxsamp(816:end,:)}; %xbackground={result_la.xxsamp(1:750,:),result_la.xxsamp(751:1500,:)};
tgrid_tc = setopt.tgrid;
xx_tc = xx;
load(['pbe_maze' num2str(maze) '_14ms.mat'])
load(['bayes_decoded_maze' num2str(maze) 'maze' num2str(maze) 'pbe.mat'])
% load('scaler_maze1_run_maze1_pbe.mat')
scaler = 0.06; %0.06; %13/500;

idx = time-time(1);
idx = round(idx./idx(2));
event_edge = event_edge+1;
%--------------
load(['sig_pbe_maze' num2str(maze) '.mat'])
pbelength=event_edge(sig_pbe_maze2_id,2)-event_edge(sig_pbe_maze2_id,1);
pbe_id = sig_pbe_maze2_id(38);
plot(xinit_ratio(event_edge(pbe_id,1):event_edge(pbe_id,2)))

n_event_range = [101,200];%[pbe_id,pbe_id];%
portion = [event_edge(n_event_range(1),1):event_edge(n_event_range(2),2)];
k=6;
% name_suffix = ['shuffle_pool_maze' num2str(12) 'pbe355136_dir.mat'];
%% ---------- run projection -----------%%
xx0=[];yy0=[];
xx0 =  [xx0;xinit_ratio(portion)];%= zeros(size(portion))'; %position(portion)';%
yy0=[yy0;double(spikes(:,portion))'];
tgrid = double(idx(portion)');%portion'; %
unique(diff(tgrid)) % pbe segment

niter = 20;
shuffle = 'no';
% figure
tic
[result, setopt, fftc, order,stepdis_ori, knndis_ori, knnportion_ori,spdcong_ori,dircong_ori] = run_projection(xx0,yy0,...
            tgrid,result_la,niter,'shuffle',shuffle,'draw',0,'xbackground',xbackground,'k',k,...
            'tctype','run','scaler',scaler,'newdatatype','pbe','givemeasure',1,'segmeasure',1,'tgrid_tc',tgrid_tc); 
        
% [result, setopt, fftc, order,] = run_projection(xx0,yy0,tgrid,result_la,niter,...
%     'shuffle',shuffle,'tctype','run','scaler',scaler,'newdatatype','pbe','draw',1); % two tracks pbe
toc

name= ['decodedPBE_maze' num2str(maze) '_20iter_' num2str(n_event_range(1)) '_' num2str(n_event_range(2)) '_' num2str(scaler) '.mat']
save(name, 'result', 'portion', 'scaler','n_event_range','setopt','xx','stepdis_ori', 'knndis_ori', 'knnportion_ori','spdcong_ori','dircong_ori');

% [result, setopt, fftc, stepdis_ori, knndis_ori, knnportion_ori,xx,order] = run_projection(xx0,yy0,tgrid,result_la,niter,'shuffle',shuffle,'xbackground',xbackground,'k',k); % two tracks
% [result, setopt, fftc, pathlen, knndis, knnportion,xx,order] =
% run_projection(xx0,yy0,tgrid,result_la,niter,'shuffle',shuffle); % run
[result, setopt, fftc] = run_projection(xx0,yy0,...
    tgrid,result_la,niter,'shuffle',shuffle,'segmeasure',0); % cross validation

%----find PBE we want ------%
d = [0,find(diff(tgrid')>1) numel(tgrid)];
i=4;

figure;
clf
% scatter(result_la.xxsamp(:,1),result_la.xxsamp(:,2),3,[.7,.7,.7])
hold on
show_latent_variable(result_la.xxsamp,xx_tc,[],tgrid_tc,'scatter_only',0)
grid on
% scatter3(result_la.xxsamp(:,1),result_la.xxsamp(:,2),result_la.xxsamp(:,3),3,[.5,.5,.5])
seg=d(i)+1:d(i+1);
% 3D
% c1=plot3(result_7ms.xxsamp(seg,1),result_7ms.xxsamp(seg,2),result_7ms.xxsamp(seg,3),'LineWidth',2,'Color',[0,1,0]);
plot3(result.xxsamp(seg,1),result.xxsamp(seg,2),result.xxsamp(seg,3),'LineWidth',2,'Color',[1,0,0]);
% plot3(result_20.xxsamp(seg,1),result_20.xxsamp(seg,2),result_20.xxsamp(seg,3),'LineWidth',2,'Color',[0,0.8,0]);

% plot(result.xxsamp(seg,1),result.xxsamp(seg,2),'LineWidth',2);
% plot(result_oldf_many.xxsamp(seg,1),result_oldf_many.xxsamp(seg,2),'LineWidth',2);
% f=find((xx_tc>xx0(seg(1))-2.5)&(xx_tc<xx0(seg(1))+2.5));
% scatter(result_la.xxsamp(f,1),result_la.xxsamp(f,2),3,'MarkerEdgeColor',[0,0,1]);
% f=find((xx_tc>xx0(seg(end))-2.5)&(xx_tc<xx0(seg(end))+2.5));
% scatter(result_la.xxsamp(f,1),result_la.xxsamp(f,2),3,'MarkerEdgeColor',[1,0,0]);

% show_latent_variable(result_la.xxsamp(seg,:),xx_tc(seg),[],tgrid_tc(seg),'scatter_only',0)
title(portion([seg(1),seg(end)]))
numel(seg)
spdcong_ori(i)
i=i+1;
%---------------------------%
figure;
show_latent_variable(setopt.xplds,xx0,xbackground,tgrid,'scatter_only',0)
title({['Maze' num2str(maze) ' first ' num2str(n_event_range(1)) '-' num2str(n_event_range(2)) ' PBEs, xinit'],['scaler=' num2str(scaler)]})
% axis([-50 50 -55 55])
% axis([-83,63,-55,50])

figure;
show_latent_variable(result.xxsamp,xx0,xbackground,tgrid,'line_only',0,'line_color',[1,0,0])
% show_latent_variable(result_la.xxsamp,xx_tc,[],tgrid_tc,'line_only',1)
% plot3(result.xxsamp(:,1),result.xxsamp(:,2),result.xxsamp(:,3),'LineWidth',2)
% show_latent_variable(result.xxsamp,xx0,[],round(tgrid),'line_only',1)
title(['Maze' num2str(maze) ' PBE No.' num2str(pbe_id) ', 20 iters'])
title({'tc: training data xx','project whole test data together','original method computing K^{1/2}'})
% axis([-50 50 -55 55])
% axis([-83,63,-55,50])


% -------- additional projection runs------- %
niter = 10;
xxtc = result_la.xxsamp;
[result,setopt] = additional_projection(xx0,yy0,fftc,xxtc,setopt,result,niter);
figure;
% xbackground={result_la.xxsamp(1:750,:),result_la.xxsamp(751:1500,:)};
show_latent_variable(result_la.xxsamp,xx_tc,[],tgrid_tc,'scatter_only',1)
show_latent_variable(result.xxsamp,xinit_ratio(portion),[],round(tgrid),'line_only',1)
title('projected x, iter = 20')


toplot = knndis_ori;
figure;
histogram(toplot(1,:),30)
hold on
histogram(toplot(2,:),30)


toplot = knnportion_ori;
figure;
histogram(toplot(1,:),30)
hold on
histogram(toplot(2,:),30)

range = 0:45;
range = 0:0.005:0.13;
figure;
histogram(toplot(1,:),range)
hold on
histogram(toplot(2,:),range)
title({['maze2 PBE No. 701-800, 14ms, 40 iters, scaler=' num2str(scaler)],'knn distance'});

%% ----------- subsequent analysis ------------%
% ---------- direction -----------%
tgrid = setopt.tgrid;
x_measure = result.xxsamp;
[speed_measure,d] = get_speed(tgrid,x_measure,1);
x_knnbase = result_la.xxsamp;
[speed_tc, d_tc] = get_speed(tgrid_tc,x_knnbase,1);
k=6;
dircong_seg = get_dir_congruence(tgrid,x_measure,x_knnbase,speed_tc,k);
[stepdis_seg, knndis_seg, knnportion_seg,spdcong_seg,dircong_seg] = projected_xx_measures_segs(x_measure,xbackground,d,tgrid,speed_tc,0,k);

% % for PBE, put in single PBE is better
% [idx_cong1,idx_cong2,congruence_ori,segcong_ori,~,~,color] = get_direction_congruence(x_measure,x_knnbase,k,speed_measure,speed_tc,d);

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

%% ---------- shuffled data pool ------------- %

% [stepdis_ori, knndis_ori, knnportion_ori] = projected_xx_measures_segs(result_30.xxsamp,xbackground,d,1,6);

pool_size = 20;
k=6;
n_manifold = numel(xbackground);
for shuftype=1:2
    if shuftype==1
        shuffle = 'cell'; %'cell';%
    else 
        shuffle = 'segtime';%
    end
    stepdis_lst = [];
    spdcong_lst = [];
    dircong_lst = [];
    for i=1:n_manifold
        knndis_lst{i}=[];
        knnportion_lst{i}=[];
    end

    tic
    for ff=1:5
        for times = 1:pool_size
            % pbe
            [~, ~, ~, ~,stepdis, knndis, knnportion,spdcong,dircong] = run_projection(xx0,yy0,...
                tgrid,result_la,niter,'shuffle',shuffle,'draw',0,'xbackground',xbackground,'k',k,...
                'tctype','run','scaler',scaler,'newdatatype','pbe','givemeasure',1,'segmeasure',1,'tgrid_tc',tgrid_tc); % two tracks pbe
            stepdis_lst = [stepdis_lst; stepdis];
            spdcong_lst = [spdcong_lst; spdcong];
            dircong_lst = [dircong_lst; dircong];
            for i=1:n_manifold
                knndis_lst{i} = [knndis_lst{i}; knndis(i,:)];
                knnportion_lst{i} = [knnportion_lst{i}; knnportion(i,:)];
            end
            % run
    %         [~, ~, ~, stepdis, knndis, knnportion,~,~] = run_projection(xx0,yy0,tgrid,result_la,niter,'shuffle',shuffle,'draw',1,'xbackground',xbackground,'k',k);
    %         stepdis_lst = [stepdis_lst, stepdis];
    %         knndis_lst = [knndis_lst, knndis];
    %         knnportion_lst = [knnportion_lst, knnportion];
        end
        save([shuffle, name_suffix],'stepdis_lst', 'knndis_lst', 'knnportion_lst','spdcong_lst','dircong_lst');
%     save([shuffle, 'shuffle_pool_maze' num2str(maze) 'pbe' num2str(pbe_id) '.mat'],'stepdis_lst', 'knndis_lst', 'knnportion_lst','spdcong_lst');
    end
    toc
end

load(['segtime',name_suffix])
timestepdis_lst=stepdis_lst;
timeknndis_lst=knndis_lst;
timeknnportion_lst=knnportion_lst;
timespdcong_lst = spdcong_lst;
timedircong_lst = dircong_lst;
load(['cell', name_suffix])

id=1;
title(['Maze' num2str(maze) ' PBE No.' num2str(n_event_range(1)-1+id) ', length: ' num2str(d(id+1)-d(id))])
measures_hist_single(id,spdcong_lst,timespdcong_lst,spdcong_ori,'speed congruence')
measures_hist_single(id,dircong_lst,timedircong_lst,dircong_ori,'direction congruence')
measures_hist_multiple(id,knnportion_lst,timeknnportion_lst,knnportion_ori,'knn portion')
measures_hist_multiple(id,knndis_lst,timeknndis_lst,knndis_ori,'knn distance')


truemaze=1;
timeknndisratio=timeknndis_lst{truemaze}./timeknndis_lst{3-truemaze};
cellknndisratio=knndis_lst{truemaze}./knndis_lst{3-truemaze};
oriknndisratio=knndis_ori(truemaze,:)./knndis_ori(3-truemaze,:);
measures_hist_single(id,cellknndisratio,timeknndisratio,oriknndisratio,['knn distance ratio, maze' num2str(truemaze) ' : maze' num2str(3-truemaze)])


figure;
hold on;
% limit = [floor(min([cellknndisratio,timeknndisratio])),ceil(max([cellknndisratio,timeknndisratio]))];
histogram(cellknndisratio,40,'Normalization','probability','BinLimits',[0,2],'DisplayStyle','stairs','EdgeColor',[0.5,0.15,.5])
histogram(timeknndisratio,40,'Normalization','probability','BinLimits',[0,2],'FaceColor',[0.5,0.15,.5],'EdgeColor','None')
plot([oriknndisratio oriknndisratio],ylim,'--','Color',[0.5,0.15,.5]*0.7,'LineWidth',2);
xlabel(['knn distance ratio, maze' num2str(truemaze) ' : maze' num2str(3-truemaze)])
ylabel('relative frequency')
legend({'cell-shuf','time-shuf','original'})


figure;scatter3(knnportion_lst{1},knnportion_lst{2},cellknndisratio,'.')
hold on
scatter3(timeknnportion_lst{1},timeknnportion_lst{2},timeknndisratio,'.')
scatter3(knnportion_ori(1),knnportion_ori(2),oriknndisratio,20,'*')
xlabel('knnportion in maze1')
ylabel('knnportion in maze2')
zlabel('knndis ratio')

figure;scatter(spdcong_lst,cellknndisratio,'.')
hold on
scatter(timespdcong_lst,timeknndisratio,'.')
scatter(spdcong_ori(1),oriknndisratio,20,'*')
xlabel('stepdis')
ylabel('knndis ratio')



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









