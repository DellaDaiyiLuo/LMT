tgrid = setopt.tgrid;
d_tc = [0 find(diff(tgrid')>1) numel(tgrid)];
xplot = result_la.xxsamp; % xinit; % xppca;%setopt.xplds;
speed = [];
for i=1:numel(d_tc)-1
    x_tmp = xplot(d_tc(i)+1:d_tc(i+1),:);
    x_sm = smoothdata(x_tmp,1,'gaussian',10);
    c = diff(x_sm);
    speed = [speed;c;c(end,:)];
end
figure;quiver(xplot(:,1),xplot(:,2),speed(:,1),speed(:,2))

xplot = [result_la.xxsamp speed];
xplot = xplot./std(xplot,[],1);

% ------
D = pdist(xplot);
Z = squareform(D);
Y = mdscale(Z,3);

%-------

k = 15;
display(['Number of neighbors to be searched: ' num2str(k)])
Idx = knnsearch(xplot,xplot,'K',k,'Distance','euclidean');

EdgeTable = zeros(size(Idx,1)*(k-1),2);
for i = 2:size(Idx,2)
    EdgeTable([size(Idx,1)*(i-2)+1:size(Idx,1)*(i-1)],:) = Idx(:,[1,i]);
end

H = simplify(graph(EdgeTable(:,1)',EdgeTable(:,2)'));
figure;plot(H)
L = laplacian(H);
[V,D] = eigs(L,5,'sa');

xx_dr2 = V(:,find(diag(D)>1e-10,3));
% figure;scatter(xx_dr2(:,1),xx_dr2(:,2),3,xx)
figure;scatter3(xx_dr2(:,1),xx_dr2(:,2),xx_dr2(:,3),3,xx)
