function logpy = estimate_py(result, fftc, xgrid, yy,tgrid,d)
nseg = numel(d)-1;
logpy = zeros(3,nseg);

for i=1:nseg
    seg = d(i)+1:d(i+1);


    ff = result.ffmat(seg,:);
    [nt,nneur] = size(ff);
    [Bfun, ~, ~, ~, ~, Kt] = prior_kernel(result.rhoxx,result.lenxx,nt,1,tgrid(seg)); 
    uu = Bfun(result.xxsamp(seg,:),1);

    logKttrm = sum(log(diag(Kt)));


    %%%%%%% cov %%%%%%%%
    covfun = covariance_fun(result.rhoff,result.lenff,2);
    cufx = covfun(xgrid,result.xxsamp(seg,:));
    cuu = covfun(xgrid,xgrid)+result.sigma*eye(size(xgrid,1));
    cuuinv = pdinv(cuu);

    % Log-determinant term
    ytrm = 0;
    etrm = 0;
    qtrm = 0;
    factrm = 0;
    for nn=1:nneur

        fn = ff(:,nn);
        yn = yy(seg,nn);
        cuuinvf_tc = cuuinv*fftc(:,nn); % (Ng,1), b
        ff_mu = cufx'*cuuinvf_tc;

        ytrm = ytrm + yn'*fn; % y*f

        etrm = etrm + sum(exp(vec(fn))); % exp(f)

        factrm = factrm + sum(log(factorial(yn)));

        fq = fn - ff_mu;

        qtrm = qtrm + fq'*fq / result.sigma;

    end

    logp_y_f = ytrm-etrm-factrm;
    logp_f_x = -.5*nt*log(result.sigma*2*pi)-.5*qtrm;
    logp_x = -.5*logKttrm-.5*nt*log(2*pi)-.5*trace(uu'*uu);

    logpy(:,i) = [logp_y_f;logp_f_x;logp_x]; 
end
end
