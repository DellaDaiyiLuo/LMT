function [result, setopt] = additional_projection(xx,yy,fftc,xxtc,setopt,result_old,niter)

nf = size(result_old.xxsamp,2);
setopt.xplds=result_old.xxsamp;
setopt.ffmat = result_old.ffmat;
setopt.niter = niter;
result = infere_latent_var(yy,nf,setopt,xx,fftc,xxtc); 
end