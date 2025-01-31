function [L,dL,ffmat] =  StateSpaceModelsofSpikeTrains_tc(ff,yymat,cufx,cuu,cuuinv,sigma2,fftc)
[nt,nneur] = size(yymat);
ffmat = reshape(ff,[],nneur);
ff = vec(ffmat);

yy = vec(yymat);
maxff = max(ff);
ff1 = ff-maxff;
log_yy_ff = yy'*ff-sum(exp(ff1))*exp(maxff);
% log_ff = -0.5*trace(ffmat'*pdinv(cufx'*cuuinv*cufx+sigma2*eye(size(cufx,2)))*ffmat);

% Quadratic term
% cuu = pdinv(cuuinv);
ff_q = ffmat - cufx'*cuuinv*fftc;
cf = cufx*ff_q;

log_ff = -.5*trace(ff_q'*ff_q)/sigma2; % posterior f cov
dL2 = -vec(ff_q/sigma2);

L = log_yy_ff+log_ff;
L = -L;

%%
dL11 = yy-exp(ff1)*exp(maxff);

dL11 = reshape(dL11,[],nneur);
dL1 = vec(dL11);
dL = dL1+dL2;
dL = -dL;

