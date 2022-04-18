% -------- load data ----------- %
% run
load('/home/kemerelab/Downloads/EEG/Elife2018_data/result_elife_t1_1000_2000_nf2_init_pca.mat')
load('/home/kemerelab/Downloads/EEG/Elife2018_data/track1_100ms.mat')
position=position(idx);
spikes=spikes(:,idx);
% simul
spikes = spikes(:,1:60)';

% -------- define input data ----------- %
tgrid_tc = setopt.tgrid;
xx_tc = xx;

portion = [2001:2785]; 

xx0 = position(portion)'; 
yy0=double(spikes(:,portion))';
tgrid = double(idx(portion)');%portion'; %

%% ---------- run projection -----------%%
niter = 20;
shuffle = 'time'; 
[result, setopt, fftc, pathlen, knndis, knnportion,xx,order] = run_projection(xx0,yy0,tgrid,result_la,niter,'shuffle',shuffle,'draw',0);

figure;
show_latent_variable(setopt.xplds,xx,{result_la.xxsamp},tgrid,[],0)
title({'Bayesian initialization of new input data'})
axis([-50 50 -55 55])

figure;
show_latent_variable(result.xxsamp,xx,{result_la.xxsamp},tgrid,[],0)
title({'New input data projected into latent space'})
axis([-50 50 -55 55])

save('4.1.2.5.mat', 'result', 'portion', 'setopt','xx','order');

% -------- additional projection runs------- %
niter = 10;
xxtc = result_la.xxsamp;
[result_20, setopt] = additional_projection(xx,yy,fftc,xxtc,setopt,result,niter);
figure;
show_latent_variable(result_20.xxsamp,xx,{result_la.xxsamp},tgrid,[],0)
title('projected x, iter = 20')


%% ----------- subsequent analysis ------------%
% ---------- direction -----------%
tgrid = setopt.tgrid;
x_measure = result.xxsamp;
[speed,d] = get_speed(tgrid,x_measure,0);
x_knnbase = result_la.xxsamp;
[speed_tc, d_tc] = get_speed(tgrid_tc,x_knnbase,0);
[congruence,segcong] = get_direction_congruence(x_measure,x_knnbase,5,speed,speed_tc,d);

figure;hold on;histogram(congruence,'FaceColor',[.5,.5,.5],'EdgeColor','None');
xlabel('speed congruence')
ylabel('count')
title('Histogram of speed congruence')
l = 400;
for i=1:numel(segcong)
    if segcong(i)>0
        p1=plot([segcong(i) segcong(i)],[0 l],'Color',[.7,.2,.2],'LineWidth',1);
    else
        p2=plot([segcong(i) segcong(i)],[0 l],'Color',[.2,.5,.7],'LineWidth',1);
    end
end
p3 = plot([0 0],[0 l],'k--','LineWidth',2);
ratio = [numel(find(segcong>0)),numel(segcong)]
if diff(ratio)==0
    legend([p1,p3],{'segment congruence>0','threshold'})
else
    legend([p1,p2,p3],{'segment congruence>0','segment congruence<=0','threshold'})
end

% ---------- shuffled data pool ------------- %
pathlen_lst = [];
knndis_lst = [];
knnportion_lst = [];

niter = 20;
shuffle = 'cell'; 
load([shuffle 'shuffle_pool.mat'])
pool_size = 200;
tic
for times = 1:pool_size
    [~, ~, ~, pathlen, knndis, knnportion,~,~] = run_projection(xx0,yy0,tgrid,result_la,niter,'shuffle',shuffle,'draw',0);
    pathlen_lst = [pathlen_lst, pathlen];
    knndis_lst = [knndis_lst, knndis];
    knnportion_lst = [knnportion_lst, knnportion];
end
toc
save([shuffle, 'shuffle_pool.mat'],'pathlen_lst', 'knndis_lst', 'knnportion_lst');
