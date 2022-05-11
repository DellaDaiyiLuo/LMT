function [result, setopt, fftc, stepdis, knndis, knnportion,xx,order] = run_projection(xx0,yy0,tgrid,result_la,niter,varargin)
% scaler: run_data * scaler = pbe_data, matrix of size (#neurons,1)
% tctype: tuning curve type, 'run' or 'pbe'
% newdatatype: new data type, 'run' or 'pbe'

% ------------ parse input ------------- %
p = inputParser;
addParameter(p,'tctype','run')
addParameter(p,'newdatatype','run')
addParameter(p,'scaler',[])
addParameter(p,'sigma',2*power(0.95,30))
addParameter(p,'draw',1)
addParameter(p,'shuffle','no')
addParameter(p,'k',10)
addParameter(p,'xbackground',[])
addParameter(p,'segmeasure',0) % get measures of every segment
parse(p,varargin{:});

scaler = p.Results.scaler;

[nt,nneur] = size(yy0);
setopt.shuffle=p.Results.shuffle;
switch p.Results.shuffle
    case 'cell'
        order = randperm(nneur);
        yy = yy0(:,order);
        xx=xx0;
    case 'time'
        order = randperm(nt);
        yy = yy0(order,:);
        xx = xx0(order);
    case 'no'
        yy=yy0;
        xx=xx0;
        order=[];
end

nf = size(result_la.xxsamp,2); % number of latent variable dimensions

if (size(tgrid,1)~=nt)||(size(yy,1)~=nt)||(size(xx,1)~=nt)
    error(['Dimensions are not cool! tgrid, xx, and yy should have ' num2str(nt) ' rows.'])
end

% --------- scaling btw tuning curve and new data --------- %

type = [p.Results.tctype p.Results.newdatatype];
switch type
    case 'runrun' % no need to scale
        fftc = result_la.ffmat;
    case 'runpbe' % scale run tc
        scaler(scaler==0) = min(nonzeros(scaler))/10;
        scale = log(scaler'); % (1,#neurons)
        fftc = result_la.ffmat+scale;
    case 'pberun' % scale run firing rate
        yy0 = yy;
        yy = zeros(size(yy0));
        for i=1:numel(scaler)
            yy(:,i) = yy0(:,i)*scaler(i);
        end
        fftc = result_la.ffmat;
    case 'pbepbe' % no need to scale
        fftc = result_la.ffmat;
end

% --------- Bayesian initialization --------- %
if nf > 2
    ng = 21;
else
    ng = 201;
end
gridbound = [];
for i=1:nf
    gridbound = [gridbound; min(result_la.xxsamp(:,i)) max(result_la.xxsamp(:,i))];
end
xgrid_ = gen_grid(gridbound,ng,nf);

fftc_init = get_tc(result_la.xxsamp,fftc,xgrid_,result_la.rhoff,result_la.lenff);

loglikelihood = -repmat(sum(exp(fftc_init),2)',nt,1) + yy*fftc_init';
[~, xinitidx] = max(loglikelihood,[],2);
xinit = xgrid_(xinitidx,:);

clear xgrid_ fftc_init gridbound loglikelihood ng xinitidx

%----use tc grid-----%
% ng = 21;
% gridbound = [];
% for i=1:nf
%     gridbound = [gridbound; min(result_la0.xxsamp(:,i)) max(result_la0.xxsamp(:,i))];
% end
% xgrid = gen_grid(gridbound,ng,nf);
% fftc = get_tc(result_la0.xxsamp,result_la0.ffmat,xgrid,result_la0.rhoff,result_la0.lenff);
%--------------------%
% run inference
setopt.xpldsmat = []; %xppcamat;%[xx xx]; %  for plotting purpose
setopt.initTYPE = 1; % initialize latent: 1. use PLDS init; 2. use random init; 3. true xx
switch setopt.initTYPE
    case 1
        setopt.xplds = xinit;% randn(nt,nf)*0.1;%%xppca;%[xx xx];% %for initialization purpose
    case 2
        setopt.xplds = [];
end

setopt.sigma2_init = p.Results.sigma; %2*power(0.95,30); % initial noise variance

setopt.rhoxx = result_la.rhoxx; % rho for Kxx
setopt.lenxx = result_la.lenxx; % len for Kxx
setopt.rhoff = result_la.rhoff; % rho for Kff
setopt.lenff = result_la.lenff; % len for Kff
setopt.la_flag = 3; % 1. no la; 2. standard la; 3. decoupled la

% same as inference
setopt.tgrid = tgrid; %[1:nt]';% ones(numel(xx)-1,1);%grid; %%double(idx)';
setopt.latentTYPE = 1; % kernel for the latent, 1. AR1 their method computing K^{-1/2}, 2. compute kernel then compute K^{-1/2}
setopt.ffTYPE = 2; % kernel for the tuning curve, 1. AR1, 2. SE
setopt.hypid = [1,2,3,4]; % 1. rho for Kxx; 2. len for Kxx; 3. rho for Kff; 4. len for Kff; 5. sigma2 (annealing it instead of optimizing it)
setopt.niter = niter; % number of iterations

% Compute P-GPLVM with Laplace Approximation
setopt.ffmat = log(yy+1e-3);
setopt.draw = p.Results.draw;
% result = infere_latent_var_xmean(yy,nf,setopt,xx,fftc,result_la.xxsamp); 
result = infere_latent_var(yy,nf,setopt,xx,fftc,result_la.xxsamp); 
d = int16([0,find(diff(tgrid')>1),numel(tgrid)]);
if p.Results.segmeasure
        [stepdis, knndis, knnportion] = projected_xx_measures_segs(result.xxsamp,p.Results.xbackground,d,p.Results.draw,p.Results.k);
else
    if isempty(p.Results.xbackground)
        [stepdis, knndis, knnportion] = projected_xx_measures(result.xxsamp,{result_la.xxsamp},d,p.Results.draw,p.Results.k);
    else
        [stepdis, knndis, knnportion] = projected_xx_measures(result.xxsamp,p.Results.xbackground,d,p.Results.draw,p.Results.k);
    end
end
end