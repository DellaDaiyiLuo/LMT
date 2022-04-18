function [result_la, setopt] = additional_pgplvm(xx,yy,setopt,result_la_old,niter)
nf = size(result_la_old.xxsamp,2);
setopt.rhoff=result_la_old.rhoff;
setopt.rhoxx=result_la_old.rhoxx;
setopt.lenff=result_la_old.lenff;
setopt.lenxx=result_la_old.lenxx;
setopt.sigma2_init=result_la_old.sigma;
setopt.xplds=result_la_old.xxsamp;
setopt.ffmat=result_la_old.ffmat;
setopt.niter = niter;
result_la = pgplvm_la(yy,nf,setopt,xx);
end