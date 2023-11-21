function L = logmargli_gplvm_se_sor_hyp_scaled_pois(loghypxx,loghypxx0,xxsamp,xgrid,latentTYPE,tgrid,nt,hypid,sigma2,ffmat,ffTYPE,cnse,yy,beta)
loghypxx0 = loghypxx;
BBwfun_x = prior_kernel(exp(loghypxx0(1)),exp(loghypxx0(2)),nt,latentTYPE,tgrid);
% imagesc(Kprior),colorbar
uu = BBwfun_x(xxsamp,1);

nneur = size(ffmat,2);
covfun = covariance_fun(exp(loghypxx0(3)),exp(loghypxx0(4)),ffTYPE); % get the covariance function
cuu = covfun(xgrid,xgrid)+cnse*eye(size(xgrid,1));

%%%%%%% cov %%%%%%%%
cufx = covfun(xgrid,xxsamp);
invcc = pdinv(cufx*cufx'+sigma2*cuu);

% Log-determinant term
logDetS1 = logdetns(cufx*cufx'+sigma2*cuu)-logdetns(cuu)+log(sigma2)*(length(xxsamp)-size(cufx,1));
logdettrm = .5*nneur*logDetS1;

% Quadratic term
cf = cufx*ffmat;
Qtrm = .5*trace(ffmat'*ffmat)/sigma2-.5*trace(invcc*cf*cf')/sigma2;

% C term: -(y*log(C)-C*exp(f)), C~N(0,beta)
C = loghypxx(5:end);
Ctrm = -sum(yy*log(C))+sum(exp(ffmat)*C); 

L = Qtrm+logdettrm+Ctrm+.5*trace(uu'*uu)+.5*beta*(C'*C);