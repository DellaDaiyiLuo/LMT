function [L,dL] = logmargli_gplvm_se_sor_la_decouple_tc(uu,yy,BBwfun,ff,covfun,sigma2,nf,BBwTfun,xgrid,cuu,cuuinv,cufx_old,invcc_old,fftc,xmean)
[nt,nneur] = size(ff);
uu = reshape(uu,[],nf);
xx = BBwfun(uu,0)+xmean; %DL

%%%%%%% cov %%%%%%%%
[cufx,dcufx] = covfun(xgrid,xx); % dcufx = d(cufx)/dx
% cuuinv = pdinv(cuu);
eff = exp(ff).*ff;

% Log-determinant term
logdetB = 0;
ytrm = 0;
etrm = 0;
qtrm = 0;
dL_logdet = 0;
dL_ytrm = 0;
dL_qtrm = 0;
dL_etrm = 0;
for nn=1:nneur
    
    fn = ff(:,nn);
    fg = fftc(:,nn);
    efn = eff(:,nn);
    yn = yy(:,nn);

    W = exp(vec(fn));
    dd = W*sigma2+1;
%     dw = sigma2+1./W;
%     invdw = 1./dw;
%     invwcf = bsxfun(@times,invdw,cufx');
%     cwc = cufx*invwcf+cuu;
%     invcwc = pdinv(cwc);
%     logdetB = logdetB+logdetns(cwc)-logdetns(cuu)+sum(log(dd));
    logdetB = logdetB+sum(log(dd)); % DL
    
    
    %%
    %     ff0 = invdd.*fn/sigma2-invdd.*cufx_old'*invcc_old*cufx_old*fn/sigma2+invdd.*efn...
    %         trm1 = -invdw*cufx'*invcwc*cufx*invdd.*fn/sigma2...
    %         trm2 = +invdw*cufx'*invcwc*cufx*invdd.*cufx_old'*invcc_old*cufx_old/sigma2*fn...
    %         trm3 = -invdw*cufx'*invcwc*cufx*invdd.*efn;
    
    invdd = 1./dd;
    
    % DL:
    cinvfg = cuuinv*fg;
    mu_old = cufx_old'*cinvfg;
    ff0 = fn - invdd.*mu_old + invdd.*(cufx'*cinvfg);
    
%     ftrm1 = cinvc(cufx,invcwc,invdw,invdd.*fn/sigma2);
%     tmp = invdd.*(cufx_old'*(invcc_old*(cufx_old*fn))/sigma2);
%     ftrm2 = cinvc(cufx,invcwc,invdw,tmp);
%     ftrm3 = cinvc(cufx,invcwc,invdw,invdd.*efn);
%     
%     ftrm0 = invdd.*fn/sigma2-tmp+invdd.*efn;
%     ff0 = ftrm0-ftrm1+ftrm2-ftrm3;
    
    %%
    % B
%     dL_logdet = dL_logdet + cinvfg*((0.5*sigma2*W.*invdd).*invdd)'; %20220224
    
    % y*f
    ytrm = ytrm+yn'*ff0;
    dL_ytrm = dL_ytrm+cinvfg*(yn.*invdd)';
    
    % exp(f)
    etrm = etrm + sum(exp(ff0));
    dL_etrm = dL_etrm+cinvfg*(exp(ff0).*invdd)';
    
    % qtrm: f'*f/K_x
    err = ff0-mu_old;
    qtrm = qtrm + err'*err/sigma2;
    dL_qtrm = dL_qtrm+cinvfg*((err/sigma2).*invdd)';
    
end

L = 0.5*logdetB-ytrm+etrm+.5*qtrm+.5*trace(uu'*uu);
% figure(3);hold on;plot(uu)


%%
dL_Kuf = dL_logdet-dL_ytrm+dL_etrm+dL_qtrm;

dL_K1 = repmat(dL_Kuf,nf,1);
dKuf1 = reshape(dcufx,[],nt);

dcf = dL_K1.*dKuf1;
dcf1 = sum(reshape(dcf,size(xgrid,1),[]),1);
dL_c = reshape(dcf1,nf,[])';

dL_u = vec(BBwTfun(dL_c,0))+vec(uu);

dL = dL_u;
