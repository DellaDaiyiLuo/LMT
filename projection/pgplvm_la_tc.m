function result = pgplvm_la_tc(yy,nf,setopt,xx,fftc,xgrid)
% Initialize the log of spike rates with the square root of spike counts.
ffmat = setopt.ffmat; %sqrt(yy);
% xllh=[];
% fllh=[];
logpy = [];


% Get sizes and spike counts
[nt,nneur] = size(yy); % nt: number of time points; nneur: number of neurons
% nf = size(xx,2); % number of latent dimensions

%
latentTYPE = setopt.latentTYPE; % kernel for the latent, 1. AR1, 2. SE
ffTYPE = setopt.ffTYPE; % kernel for the tuning curve, 1. AR1, 2. SE
if nf==1
    xpldsmat = setopt.xpldsmat;
end
xplds = setopt.xplds;

% generate grid values as inducing points
tgrid = setopt.tgrid;

% set hypers
hypers = [setopt.rhoxx, setopt.lenxx, setopt.rhoff, setopt.lenff]; % rho for Kxx; len for Kxx; rho for Kff; len for Kff

% set initial noise variance for simulated annealing
% lr = setopt.lr; % learning rate
sigma2 = setopt.sigma2_init;
propnoise_init = 0.001;
propnoise = propnoise_init;

% set initial prior kernel
% K = Bfun(eye(nt),0)*Bfun(eye(nt),0)';
% Bfun maps the white noise space to xx space
[Bfun, BTfun, nu, sdiag, iikeep, Kprior] = prior_kernel(hypers(1),hypers(2),nt,latentTYPE,tgrid);
rhoxx = hypers(1); % marginal variance of the covariance function the latent xx
lenxx = hypers(2); % length scale of the covariance function for the latent xx
rhoff = hypers(3); % marginal variance of the covariance function for the tuning curve ff
lenff = hypers(4); % length scale of the covariance function for the tuning curve ff


result.rhoxx = rhoxx;
result.lenxx = lenxx;
result.rhoff = rhoff;
result.lenff = lenff;
result.sigma = sigma2;

% initialize latent
initTYPE = setopt.initTYPE;
switch initTYPE
    case 1  % use LLE or PPCA or PLDS init
        uu0 = Bfun(xplds,1);
    case 2   % use random init
        uu0 = randn(nu,nf)*0.01;
    case 3   % true xx
        uu0 = Bfun(xx,1);
end
uu = uu0;  % initialize sample
xxsamp = Bfun(uu,0);
if nf==1
    xxsampmat = align_xtrue(xxsamp,xpldsmat);
    xxsampmat_old = xxsampmat;
    xpldsmat = xxsampmat;
end

% Now do inference
infTYPE = 1; % 1 for MAP; 2 for MH sampling; 3 for hmc
ppTYPE = 1; % 1 optimization for ff; 2. sampling for ff
la_flag = setopt.la_flag; % 1. no la; 2. standard la; 3. decoupled la
% opthyp_flag = setopt.opthyp_flag; % flag for optimizing the hyperparameters

% set options for minfunc
options = [];
options.Method='scg';
options.TolFun=1e-4;
options.MaxIter = 1e1;
options.maxFunEvals = 1e1;
options.Display = 'off';

niter = setopt.niter;
clf
for iter = 1:niter
%     display(['iter' num2str(iter)])
    
    %% 1. Find optimal ff
     covfun = covariance_fun(rhoff,lenff,ffTYPE); % get the covariance function
     cuu = covfun(xgrid,xgrid)+sigma2*eye(size(xgrid,1));
     cuuinv = pdinv(cuu);
     cufx = covfun(xgrid,xxsamp);
     
     lmlifun_poiss = @(ff) StateSpaceModelsofSpikeTrains_tc(ff,yy,cufx,cuuinv,sigma2,fftc);

    
    switch ppTYPE
        case 1
            ff0 = vec(ffmat);
            floss_ff = @(ff) lmlifun_poiss(ff); % negative marginal likelihood
            % DerivCheck(floss_ff,ff0)
            [ffnew, fval] = minFunc(floss_ff,ff0,options);
        case 2
            % set up MCMC inference
            nsperiter = 10;
            fproprnd_ff = @(ff)(ff+randn(size(ff))*0.1); % proposal distribution
            flogpdf_ff = @(ff)(-lmlifun_poiss(ff'));
            
            ff0 = vec(Bfun_cov(ffmat,1));
            ffnew = mhsample_anqi(ff0',nsperiter,'logpdf',flogpdf_ff,'proprnd',fproprnd_ff,'symmetric',true)';
            ffnew = ffnew(:,end);
    end
    [L,dL,ffnew] = lmlifun_poiss(ffnew);
%     display(['optimal ff loss: ' num2str(L)])
%     fllh = [fllh,L];
    ffmat = ffnew;
    
%     subplot(412),plot([ff(:,20),exp(ffnew(:,20))]),title('exp(ff)'),legend('true ff','P-GPLVM ff'),drawnow
    
    %% 2. Find optimal latent xx, actually search in u space, xx=K^{1/2}*u
    [Bfun, BTfun, nu] = prior_kernel(rhoxx,lenxx,nt,latentTYPE,tgrid); % DL: using Cholesky decomposition to compute K^{1/2} and K^{1/2}.T, getting uu~N(0,1)
    uu = Bfun(xxsamp,1);
    cufx_old = covfun(xgrid,xxsamp);
    invcc_old = pdinv(cufx_old*cufx_old'+sigma2*cuu);
    
    switch ffTYPE
        case 1 % AR1 without grad
            % lmlifun = @(u) logmargli_gplvm_ar(u,Bfun,ffmat,covfun,sigma2,nf); % only works for 1d
            lmlifun = @(u) logmargli_gplvm_se(u,Bfun,ffmat,covfun,sigma2,nf);
        case 2 % SE with grad
            switch la_flag
                case 1
                    % no la, poisson
                    lmlifun = @(u) logmargli_gplvm_se_sor(u,Bfun,ffmat,covfun,sigma2,nf,BTfun,xgrid,cuu);
                case 2
                    % standard la
                    lmlifun = @(u) logmargli_gplvm_se_sor_la(u,Bfun,ffmat,covfun,sigma2,nf,BTfun,xgrid,cuu);
                case 3
                    % decouple la
                    lmlifun = @(u) logmargli_gplvm_se_sor_la_decouple_tc(u,yy,Bfun,ffmat,covfun,sigma2,nf,BTfun,xgrid,cuu,cuuinv,cufx_old,fftc);
            end
    end
    
    % set up MAP inference
    floss = @(u) lmlifun(vec(u));
    opts = optimset('largescale', 'off', 'maxiter', 15, 'display', 'iter');
    
    % set up MCMC inference
    nsperiter = 50;
    fproprnd = @(u)(u+randn(size(u))*propnoise); % proposal distribution
    flogpdf = @(u) -lmlifun(u');
    
    % set up HMC inference, if use HMC, we have to choose the SE kernel
    % which returns grad
    flogpdf_grad = @(u) hmc_grad(u,lmlifun);
    
    % ========================================
    switch infTYPE
        case 1, % do MAP infernece
            switch ffTYPE
                case 1 % AR1, fminunc, no grad
                    uunew = fminunc(floss,vec(uu),opts);
                case 2 % SE, minFunc, with grad
                    % DerivCheck(floss,vec(randn(size(uu))))
                    [uunew, fval] = minFunc(floss,vec(uu),options);
            end
        case 2, % do MCMC
            uunew = mhsample(vec(uu)',nsperiter,'logpdf',flogpdf,'proprnd',fproprnd,'symmetric',true)';
            uunew = uunew(:,end);
        case 3, % do HMC
            options_hmc = foptions;             % Default options vector.
            options_hmc(1) = 1;         % Switch on diagnostics.
            options_hmc(7) = 1;     % Number of steps in trajectory.
            options_hmc(14) = 1;        % Number of Monte Carlo samples returned.
            options_hmc(15) = 10;       % Number of samples omitted at start of chain.
            options_hmc(18) = 0.001;        % Step size.
            
            % hmc('state', 42);
            uunew = hmc(floss, vec(uu)', options_hmc, flogpdf_grad)';
            uunew = uunew(:,end);
    end
    uu = reshape(uunew,[],nf);
    xxsamp = Bfun(uu,0);
%     display(['optimal xx loss: ' num2str(fval)])
%     xllh=[xllh,fval];
    
    % plot latent xx
    hold on
        if setopt.draw
        switch nf
            case 1
                xxsampmat = align_xtrue(xxsamp,xx);
                subplot(212); plot(1:nt,xx,'b-',1:nt,xpldsmat,'m.-',1:nt,xxsampmat,'k-',1:nt,xxsampmat_old,'k:','linewidth',2); legend('true x','init x','P-GPLVM x','P-GPLVM old x');
                xlabel('time bin'); drawnow;
                xxsampmat_old = xxsampmat;
            case 2
                scatter(xxsamp(:,1),xxsamp(:,2),3,'MarkerEdgeColor',[1/niter*iter,0,0]); drawnow;
            case 3
                show_latent_variable(xxsamp,xx,[],tgrid,'line_only',1,'line_color',[1/niter*iter,0,0])
%                 scatter3(xxsamp(:,1),xxsamp(:,2),xxsamp(:,3),3,'MarkerEdgeColor',[1/niter*iter,0,0]); 
                drawnow;
        end
        end

%     logpy = [logpy,sum(estimate_py(result, fftc, xgrid, yy,tgrid,d))];
end
% figure;plot(xllh);title('x loss')
% figure;plot(fllh);title('f loss')
% figure;plot(logpy);title('logpy')

result.ffmat = ffmat;
result.xxsamp = xxsamp;
% result.xloss = xllh(end);
% result.floss = fllh(end);



