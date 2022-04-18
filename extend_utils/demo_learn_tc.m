
initpaths;
addpath(genpath('/home/kemerelab/Downloads/EEG/utils'));

%% Load data

datatype = 'simul';

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
        plot(position(idx));
        position = position(idx);
        spikes = spikes(:,idx);
    case 'pbe'
        load('track1_PBEs_4ms.mat')
        load('pbe_bayesian_decoded_ratio_runall.mat') 
        position = xinit_ratio'; % (optional) just for coloring the figure
        clear xinit_ratio
end

portion = [2001:3500];% 671:820];% 901:1050];
xx = position(portion)'; 
figure;scatter(1:numel(portion),xx,3,xx)
nf = 2;  % number of latent dimensions
niter = 30;
yy = double(spikes(:,portion))';
tgrid = portion';%double(idx(portion))';%[1:numel(portion)]'; %

[result_la,setopt]=run_pgplvm(xx,yy,tgrid,nf,niter,[]);
figure
show_latent_variable(result_la.xxsamp,xx,[],tgrid,[],1)

% additional run if not satisfied
niter=20;
[result_la,setopt] = additional_pgplvm(xx,yy,setopt,result_la,niter);

save('result_4.1.2_30iter.mat', 'result_la', 'xx', 'setopt','portion');


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
direct=zeros(nt,1);direct(diff(xx)<0)=1;
direct(end) = direct(end-1);

    % show directions in position
figure;hold on;
scatter(ii(direct==1),xx(ii(direct==1)),3);
scatter(ii(direct==0),xx(ii(direct==0)),3)

    % show directions in latent space
figure;hold on;scatter(xplot(direct==1,1),xplot(direct==1,2),3)
scatter(xplot(direct==0,1),xplot(direct==0,2),3)
legend({'Forward','Backward'})
axis([-50,50,-55,55])

% ---------- number of manifolds ---------- %
ratio = 0.002;
xplot = result_la.xxsamp;%setopt.xplds;%

[tbl_id, comp_id,k_log,n_debris_log] = find_components(xplot, ratio, 1);


figure;hold on
n_comp = size(tbl_id,1);
title_cell = cell(n_comp);
for i=1:n_comp
    show_latent_variable(xplot,comp_id,[],[],comp_id==tbl_id(i,1),1)
    title_cell(i) = ['comp' num2str(i)];
end

legend(title_cell)



