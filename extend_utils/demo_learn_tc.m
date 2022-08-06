
addpath(genpath('/home/kemerelab/Downloads/EEG/LMT-master'));warning off
addpath(genpath('/home/kemerelab/Downloads/EEG/extend_utils'));

%% Load data

datatype = 'run';
env = 2;
if env<2
    switch datatype % spikes: (#neuron,#sample)
        case 'simul'
            load('simulated3_spike_track1_60dir.mat')
            spikes0 = spikes;
            cell_id = [1:60];
            pos = linspace(0,10,201);
            figure;hold on;
            for i=cell_id
                p2=plot(pos,p_bwd(:,i)+i,'k');
                p1=plot(pos,p_fwd(:,i)+i,'Color',[.7,0,0]);
            end
            xlim([0,10])
            legend([p1,p2],{'Forwad','Backward'})
            spikes = spikes0(:,cell_id)';
        case 'run'
            load('track1_100ms.mat')
            idx=idx+1;
            plot(pos_linear(idx));
            position = pos_linear(idx);
            spikes = spikes(:,idx);
        case 'pbe'
            load('track1_PBEs_4ms.mat')
            load('pbe_bayesian_decoded_ratio_runall.mat') 
            position = xinit_ratio'; % (optional) just for coloring the figure
            clear xinit_ratio
    end
    portion = [1:891];% 671:820];% 901:1050];
    xx = position(portion); 
    figure;scatter(1:numel(portion),xx,3,xx)
    yy = double(spikes(:,portion))';
    tgrid = double(idx(portion))';%[1:numel(portion)]'; %portion';%
else
    load('maze1_run_v1_3_v2_5_500ms_new.mat')
    position1 = pos_linear;
    plot(position1);
    idx = time-time(1);
    idx = round(idx./idx(2));
    spikes1 = double(spikes);
    tgrid1 = double(idx)';
    load('maze2_run_v1_3_v2_5_500ms_new.mat')
    position2 = pos_linear;
    plot(position2);
    idx = time-time(1);
    idx = round(idx./idx(2));
    spikes2 = double(spikes);
    tgrid2 = double(idx)';
    
    xx = [position1,position2]';
    yy = [spikes1,spikes2]';
    tgrid = [tgrid1;tgrid2+max(tgrid1)+1000];
    plot(xx);
end


figure
nf = 3;  % number of latent dimensions
niter = 30;
[result_la,setopt]=run_pgplvm(xx,yy,tgrid,nf,niter,[]);
show_latent_variable(result_la.xxsamp,xx,[],setopt.tgrid,'scatter_only',0)

% additional run if not satisfied
niter=20;
[result_la,setopt] = additional_pgplvm(xx,yy,setopt,result_la,niter);

% save('result_la_maze2_3nf_30iter_idxtgrid.mat', 'result_la', 'xx', 'setopt','portion1','portion2');
save('result_la_maze12_500ms_3nf_30iter_allrun_idxtgrid.mat', 'result_la', 'xx', 'setopt');



%% After inference
%% simulation
% --------- compare tuning curves ----------- %
ng = 201;
xgrid_ = gen_grid([min(result_la.xxsamp) max(result_la.xxsamp)],ng,nf);

    % nondirect 1D tc
fftc = get_tc(result_la.xxsamp,result_la.ffmat,xgrid_,result_la.rhoff,result_la.lenff);
fftc = fftc(end:-1:1,:);
figure;hold on;
for i=1:40
    m = max(p(:,i))/2;
    plot(p(:,i)/m+i,'k')
    plot(exp(fftc(:,i))/m+i,'r')
end
axis([0,200,0,43])

%% run
% ---------- directions ---------- %
xplot = result_la.xxsamp;
nt=size(xx,1);
ii=1:nt;
direct=zeros(nt-1,1);direct(diff(xx)<0)=1;
direct = [direct(1);direct];
direct_sm = smoothdata(direct,'movmedian',5);
direct = round(direct_sm);

    % show directions in position
figure;hold on;
scatter(ii(direct==1),xx(ii(direct==1)),3);
scatter(ii(direct==0),xx(ii(direct==0)),3)
    % show directions in latent space
figure;hold on;scatter(xplot(direct==1,1),xplot(direct==1,2),3)
scatter(xplot(direct==0,1),xplot(direct==0,2),3)
legend({'Forward','Backward'})
axis([-83,63,-55,50])
axis([-50,50,-55,55])

% ---------- number of manifolds ---------- %
ratio = 0.002;
xplot = result_la.xxsamp;%([1:815,1216:end],:);%setopt.xplds;%

[tbl_id, comp_id,k_log,n_debris_log] = find_components(xplot, ratio, 1);
figure;plot(k_log,n_debris_log)

n_comp = size(tbl_id,1);
xbackground = cell(n_comp,1);
for i=1:n_comp
    xbackground{i} = result_la.xxsamp(comp_id==tbl_id(i,1),:);
end
k = k_log(end);
save('/home/kemerelab/Downloads/EEG/Elife2018_data/True_Shuffle_MIX/result_elife_mix_nf2_init_pca_xbackground.mat', 'xbackground', 'k')

figure
ax1 = axes;
i=1;
scatter(ax1,xbackground{i}(:,1),xbackground{i}(:,2),3,xx(comp_id==tbl_id(i,1)))
view(2)
ax2 = axes;
i=2;
scatter(ax2,xbackground{i}(:,1),xbackground{i}(:,2),3,xx(comp_id==tbl_id(i,1)))
%%Link them together
linkaxes([ax1,ax2])
%%Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'parula')
colormap(ax2,'copper')
axis([-83,63,-55,50])
%%Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 0.036 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 0.036 .815]);




