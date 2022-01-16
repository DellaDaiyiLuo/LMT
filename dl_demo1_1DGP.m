% cd /home/kemerelab/Downloads/EEG/LMT-master
% demo1_1DGP.m
%
% Tutorial script illustrating P-GPLVM for 1-dimensional latent variable
% with tuning curves generated from 1D Gaussian Process.

% Initialize paths
initpaths;

% Load data
load('/home/kemerelab/Downloads/EEG/simulated3_spike_maze2.mat')
load('maze2_run.mat')
load('simulated3_spike_maze2_60.mat', 'fwd')
load('simulated3_spike_maze2_60.mat', 'bwd')

figure;plot(pos_linear);hold on;plot(fwd,pos_linear(fwd),'.');plot(bwd,pos_linear(bwd),'.');


fwd = fwd+1;
bwd = bwd+1;
fwdidx = fwd(fwd>=1000 & fwd<=2000)-999;
bwdidx = bwd(bwd>=1000 & bwd<=2000)-999;

pos = (pos-min(pos))/(max(pos)-min(pos));
pos = pos * 10;

xx = position(1000:2000)';
yy = double(spikes(1000:2000,1:40));
% yy = double(spikes(:,1000:2000))';

% plot tuning curve
figure;hold on; 
k=0;
xgrid = linspace(0,10,201);
for i = 41:60
    plot(xgrid,p(:,i)/2+k,'r')
%     plot(ppp,p_bwd(:,i)/2+k,'k')
    k=k+1;
end
% ff = simdata.spikeRates;

% Get sizes and spike counts
[nt,nneur] = size(yy); % nt: number of time points; nneur: number of neurons
nf = 1;  % number of latent dimensions

setopt.tgrid = [1:nt]';% ones(numel(xx)-1,1);%grid; %double(idx(1000:2000)'); %double(idx(1000:2000))';%
setopt.latentTYPE = 1; % kernel for the latent, 1. AR1 their method computing K^{-1/2}, 2. compute kernel then compute K^{-1/2}

% plot(position);hold on;plot(bwd,xx,'.')
%% == 1. Compute baseline estimates ====

% Initialize the log of spike rates with the square root of spike counts.
ffmat = sqrt(yy);

% % Compute LLE
% xlle = lle(ffmat,nf,9);
% xllemat = align_xtrue(xlle,xx); % align the estimate with the true latent variable.
% plot(xx);hold on;plot(xllemat); % scatter(xlle(:,1),xlle(:,2),3,xx)

% % Compute PPCA
xppca = pca(ffmat,nf);
% xppca = [xppca;pca(ffmat,nf)];
xppcamat = align_xtrue(xppca,xx);
% options = fgplvmOptions('ftc');
% xppca = genX_ppca(nf, nneur, ffmat, options);
% xppcamat = align_xtrue(xppca,xx); % align the estimate with the true latent variable.

% Compute Poisson Linear Dynamic System (PLDS)
% xplds = run_plds(yy,nf)';
% xpldsmat = align_xtrue(xplds,xx); % align the estimate with the true latent variable.
% xpldsmat = xppcamat;
setopt.xpldsmat = xppcamat;%[xx xx]; %xx; %  for plotting purpose

% xplds = xppca;
% plot(xx);hold on;plot(xpldsmat); % 
% scatter(xplds(:,1),xplds(:,2),3,xx)

% grid = double(diff(idx))';
% tgrid = ones(size(grid));
% tgrid(grid>200)=inf;

%% == 2. Compute P-GPLVM ====
% Set up options

setopt.lr = 0.95; % learning rate
setopt.ffTYPE = 2; % kernel for the tuning curve, 1. AR1, 2. SE

setopt.initTYPE = 1; % initialize latent: 1. use PLDS init; 2. use random init; 3. true xx
setopt.opthyp_flag = 1;
setopt.sigma2_init = 2; % initial noise variance
setopt.rhoxx = 10; % rho for Kxx
setopt.lenxx = 100; % len for Kxx
setopt.rhoff = 10; % rho for Kff
setopt.lenff = 50; % len for Kff


setopt.la_flag = 3; % 1. no la; 2. standard la; 3. decoupled la
setopt.hypid = [1,2,3,4]; % 1. rho for Kxx; 2. len for Kxx; 3. rho for Kff; 4. len for Kff; 5. sigma2 (annealing it instead of optimizing it)
%setopt.xpldsmat = xllemat; % for plotting purpose
%setopt.xplds = xlle; % for initialization purpose
setopt.xplds = xppca;%[xx xx];%xx; % %for initialization purpose
setopt.niter = 20; % number of iterations



% Compute P-GPLVM with Laplace Approximation
result_la = pgplvm_la(yy,nf,setopt,xx);

% additional run
setopt.rhoff=result_la.rhoff;
setopt.rhoxx=result_la.rhoxx;
setopt.lenff=result_la.lenff;
setopt.lenxx=result_la.lenxx;
setopt.sigma2_init=result_la.sigma;
setopt.xplds=result_la.xxsamp;
setopt.niter = 10;
result_la = pgplvm_la(yy,nf,setopt,xx);


ii=1:1001;
direct=zeros(1001,1);direct(diff(xx)<0)=1;

figure;plot(xx);hold on;plot(ii(direct==1),xx(direct==1));plot(ii(direct==0),xx(direct==0))
% xxsampmat = align_xtrue(result_la.xxsamp,xx);
% figure;plot(align_xtrue(result_la.xxsamp(:,1),xx));hold on;plot(align_xtrue(result_la.xxsamp(:,2),xx));plot(xx)
% 
% xgrid = gen_grid([min(result_la.xxsamp(:,1)) max(result_la.xxsamp(:,1))],201,nf); % x grid for plotting tuning curves
% fftc = exp(get_tc(result_la.xxsamp,result_la.ffmat,xgrid,result_la.rhoff,result_la.lenff));

%3d
figure;scatter3(result_la.xxsamp(:,1),result_la.xxsamp(:,2),result_la.xxsamp(:,3),3,xx)

figure;scatter3(result_la.xxsamp(fwdidx,1),result_la.xxsamp(fwdidx,2),result_la.xxsamp(fwdidx,3),10,'filled')
hold on;
scatter3(result_la.xxsamp(bwdidx,1),result_la.xxsamp(bwdidx,2),result_la.xxsamp(bwdidx,3),10,'filled')


%2d
figure;scatter(result_la.xxsamp(:,1),result_la.xxsamp(:,2),3,xx)

figure;scatter(result_la.xxsamp(fwdidx,1),result_la.xxsamp(fwdidx,2),10,'filled')
hold on;
scatter(result_la.xxsamp(bwdidx,1),result_la.xxsamp(bwdidx,2),10,'filled')


save('result_sim3_spike_maze2_nond_60_cell_1_40.mat', 'result_la', 'xx', 'setopt');

% Compute P-GPLVM with a variational lower bound
% result_va = pgplvm_va(yy,xx,setopt);

%% == 3. Plot latent variables and tuning curves ====

% 1D tuning curve
xxsampmat = align_xtrue(result_la.xxsamp,xx);
xgrid = gen_grid([min(xxsampmat(:,1)) max(xxsampmat(:,1))],201,nf); % x grid for plotting tuning curves
fftc = get_tc(xxsampmat,result_la.ffmat,xgrid,result_la.rhoff,result_la.lenff);
yytc = exp(fftc);

% 2D to 1D tuning curve
tbl = tabulate(position);
pp = round(position*20);
xx_ = (xx-min(xx))/(max(xx)-min(xx))*10;
portion = direct==0; %bwdidx; % fwd/bwd/all
x = round(xx_(portion)*20);
xxsamp = result_la.xxsamp(portion,:);
tbl1 = tabulate(x);
xgrid = zeros(size(tbl1,1), nf);
xgrid_ = zeros(201, nf);
for i=1:size(tbl1,1)
    xgrid(i,:)=mean(xxsamp(x==tbl1(i,1),:),1);
end
for i=1:nf
    xgrid_(:,i)=smooth(interp1(tbl1(:,1),xgrid(:,i),0:200,'nearest'),10);
end

xgrid_1=xgrid_;
xgrid_0=xgrid_;

figure;
if nf==2
    % dir
    scatter(result_la.xxsamp(direct==0,1),result_la.xxsamp(direct==0,2),3);
    hold on;
    plot(xgrid_0(:,1),xgrid_0(:,2),'k-')
    
    % nondir
%     scatter(result_la.xxsamp(:,1),result_la.xxsamp(:,2),3,xx);
%     hold on;
%     plot(xgrid_(:,1),xgrid_(:,2),'r-')
elseif nf==3
    scatter3(result_la.xxsamp(:,1),result_la.xxsamp(:,2),result_la.xxsamp(:,3),3,xx)
    hold on;
    plot3(xgrid_(:,1),xgrid_(:,2),xgrid_(:,3),'k-')
end

% infer tuning curve
fftc = get_tc(result_la.xxsamp(portion,:),result_la.ffmat(portion,:),xgrid_,result_la.rhoff,result_la.lenff);
yytc = exp(fftc);

fftc1 = fftc;
fftc0 = fftc;



% plot tunning curve
figure;
% yytc = exp(fftc);
yytc1 = exp(fftc1);
yytc0 = exp(fftc0);
% p = p_bwd(:,cell_id);
for i=1:40
    subplot(6,8,i)
    hold on;
    plot(linspace(0,10,201),p(:,i));
    plot(linspace(0,10,201),yytc(:,i));title(['Cell', num2str(i)]);
%     plot(linspace(0,10,201),yytc1(:,i));title(['Cell', num2str(i)]);
%     plot(linspace(0,10,201),yytc0(:,i));title(['Cell', num2str(i)]);
end
legend('dir1','dir0')

save('result_simul_maze2_1000_2000_nf1_40_init_pca.mat', 'result_la', 'xx', 'setopt','fftc','xgrid');
% %%
% 
% plot(xx);hold on;plot(xxsampmat);
% figure
% scatter(xx, xxsampmat);
% scatter(result_la.xxsamp(:,1),result_la.xxsamp(:,2),3,xx)
% 
% 
% subplot(211); plot(1:nt,xx','b-',1:nt,xpldsmat,'m.-',1:nt,xxsampmat,'k:','linewidth',2); legend('true x','LLE x','PLDS x','P-GPLVM x');
% xlabel('time bin'); drawnow;

% gridbound = [];
%     for i=1:nf
%         gridbound = [gridbound; min(result_la.xxsamp(:,i)) max(result_la.xxsamp(:,i))];
%     end
% xgrid = gen_grid(gridbound,25,nf); % x grid for plotting tuning curves
% fftc = exp(get_tc(result_la.xxsamp,result_la.ffmat,xgrid_,result_la.rhoff,result_la.lenff));

% neuronlist = 1:nneur;
% ii = randperm(nneur);
% neuronlist = neuronlist(ii(1:4));
% for ii=1:4
%     neuronid = neuronlist(ii);
%     subplot(4,2,4+ii), cla
%     hold on,
%     plot(simdata.xgrid,simdata.tuningCurve(:,neuronid),'b-')
%     plot(xgrid,fftc(:,neuronid),'r-')
%     legend('true tc','estimated tc')
%     title(['neuron ' num2str(neuronid)])
%     hold off
% end
% 


%% predict position

load('result_simul_sin2_20_nf1_init_true.mat')
load('simulated3_spike_sine_20.mat')
yy = double(spikes(:,:))';
xx = position(:); 
[nt,nneur] = size(yy); % nt: number of time points; nneur: number of neurons


nf=1;
% tc for initialization
ng = 201;
gridbound = [];
for i=1:nf
    gridbound = [gridbound; min(result_la.xxsamp(:,i)) max(result_la.xxsamp(:,i))];
end
xgrid_ = gen_grid(gridbound,ng,nf);
% nondirect 1D tc
if nf==1
    xxsampmat = align_xtrue(result_la.xxsamp,xx);
    fftc_init = get_tc(xxsampmat,result_la.ffmat,xgrid_,result_la.rhoff,result_la.lenff);
else
    fftc_init = get_tc(result_la.xxsamp,result_la.ffmat,xgrid_,result_la.rhoff,result_la.lenff);
end
% xgrid_=xgrid;
loglikelihood = -repmat(sum(exp(fftc_init),2)',nt,1) + yy*fftc_init';
[~, xinitidx] = max(loglikelihood,[],2);
xinit = xgrid_(xinitidx,:);
switch nf
    case 1
        plot(xx);hold on;plot(xinit);
    case 2
        scatter(xinit(:,1),xinit(:,2),3,xx);
end

% direct tc
d = 1;
yy_dir = yy(direct==d,:);
xx_dir = xx(direct==d);
fftc_init = get_tc(result_la.xxsamp(direct==d,:),result_la.ffmat(direct==d,:),xgrid_,result_la.rhoff,result_la.lenff);
loglikelihood = -repmat(sum(exp(fftc_init),2)',size(xx_dir,1),1) + yy_dir*fftc_init';
[~, xinitidx] = max(loglikelihood,[],2);
xinit = xgrid_(xinitidx,:);
plot(xx_dir);hold on;plot(xinit)

% plot xinit
figure;
switch nf
    case 1
        plot(position);hold on;plot(xinit);
    case 2
        scatter(xinit(:,1),xinit(:,2),3,xx_dir);
end

% tc for running inference
switch nf
    case 1
        ng = 201;
    case 2
        ng = 25;
    case 3
        ng = 10;
end
gridbound = [];
for i=1:nf
    gridbound = [gridbound; min(result_la.xxsamp(:,i)) max(result_la.xxsamp(:,i))];
end
xgrid = gen_grid(gridbound,ng,nf);
fftc = get_tc(result_la.xxsamp,result_la.ffmat,xgrid,result_la.rhoff,result_la.lenff);
fftc1 = get_tc(result_la.xxsamp(direct==1,:),result_la.ffmat(direct==1,:),xgrid,result_la.rhoff,result_la.lenff);

fftc0 = get_tc(result_la.xxsamp(direct==0,:),result_la.ffmat(direct==0,:),xgrid,result_la.rhoff,result_la.lenff);

% different from inference
setopt.xpldsmat = xinit; %xppcamat;%[xx xx]; %  for plotting purpose
setopt.xplds = xinit; %xppca;%[xx xx];% %for initialization purpose
setopt.opthyp_flag = 0;
setopt.sigma2_init = 2; %result_la.sigma; % initial noise variance
setopt.rhoxx = result_la.rhoxx; % rho for Kxx
setopt.lenxx = result_la.lenxx; % len for Kxx
setopt.rhoff = result_la.rhoff; % rho for Kff
setopt.lenff = result_la.lenff; % len for Kff
setopt.la_flag = 3; % 1. no la; 2. standard la; 3. decoupled la

setopt.initTYPE = 1; % initialize latent: 1. use PLDS init; 2. use random init; 3. true xx


% same as inference
setopt.tgrid = double(idx'); %[1:nt]';% ones(numel(xx)-1,1);%grid; %%double(idx)';
setopt.latentTYPE = 1; % kernel for the latent, 1. AR1 their method computing K^{-1/2}, 2. compute kernel then compute K^{-1/2}
setopt.lr = 1; % learning rate
setopt.ffTYPE = 2; % kernel for the tuning curve, 1. AR1, 2. SE
setopt.hypid = [1,2,3,4]; % 1. rho for Kxx; 2. len for Kxx; 3. rho for Kff; 4. len for Kff; 5. sigma2 (annealing it instead of optimizing it)
setopt.niter = 15; % number of iterations

nf=1;

% Compute P-GPLVM with Laplace Approximation
result = infere_latent_var(yy,nf,setopt,xx,fftc,xgrid);
