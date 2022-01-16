function L = logmargli_gplvm_se_sor_hyp(loghypxx,loghypxx0,seg,yymat,xgrid,latentTYPE,nt,hypid,sigma2,fgrid,ffTYPE,ntr)
loghypxx0(hypid) = loghypxx;

%%
nf = size(seg(1).xxsamp,2);
uu = [];
for i=1:ntr
    [Bfun, BTfun, nu] = prior_kernel_sp(exp(loghypxx0(1)),exp(loghypxx0(2)),nt(i),latentTYPE);
    uu = [uu; Bfun(seg(i).xxsamp,1)];
end


%%%%%%% cov %%%%%%%%
covfun = covariance_fun(exp(loghypxx0(3)),exp(loghypxx0(4)),ffTYPE); % get the covariance function
[kxx,dcc] = covfun(xgrid,xgrid);
sigma2 = kxx(1,1)*sigma2;
invkxx = pdinv(kxx+sigma2*eye(size(kxx)));

% poisson
xxsamp_mt = [];
for i=1:ntr
    xxsamp_mt = [xxsamp_mt; seg(i).xxsamp];
end
ctx = covfun(xgrid,xxsamp_mt);
ctx = ctx';
invkf = invkxx*fgrid;
ffmat = ctx*invkf;

ff = vec(ffmat);
yy = vec(yymat);
log_yy_ff = yy'*ff-sum(exp(ff));

L = -log_yy_ff+.5*trace(uu'*uu);
