function [result, setopt] = additional_projection(xx,yy,fftc,xxtc,setopt,result_old,niter)

nf = size(result_old.xxsamp,2);
setopt.xplds=result_old.xxsamp;
setopt.ffmat = result_old.ffmat;
setopt.niter = niter;
result = pgplvm_la_tc(yy,nf,setopt,xx,fftc,xxtc);
end