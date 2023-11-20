function show_latent_variable(xplot,xcolor,xbackground,tgrid,varargin)
% nf: scalar, latent variable dimensions
% xplot: (nt,nf) matrix, latent variable to be plotted
% xcolor: (nt,1) matrix, color indices for xplot, usually true animal position
% xbackground: *cells* of (n,nf) matrices, background latent variables
% (showing tuning grid)
% tgrid: showing continuouty of xplot
% part: showing part of xplot
% scatter_only: only scatter xplot, no segment connection shown
input = inputParser;
addParameter(input,'part',1:size(xplot,1))
addParameter(input,'scatter_only',0)
addParameter(input,'line_only',0)
addParameter(input,'line_color',[0.7,0.5,0.5])
addParameter(input,'legend',1)
parse(input,varargin{:});

ii = 1:size(xplot,1);
part = ii(input.Results.part);

nf = size(xplot,2);
hold on

switch nf
    case 1
        xxsampmat = align_xtrue(xplot,xcolor);
        plot(xcolor)
        plot(xxsampmat)
        legend({'true xx','estimated xx'})
    case 2
        p = [];
        label = {};
        if ~isempty(xbackground)
            n_xbg = numel(xbackground);
            for i=1:n_xbg
                p1=scatter(xbackground{i}(:,1),xbackground{i}(:,2),50,'MarkerEdgeColor',[.8,.8,.8]*i/n_xbg,'Marker','.');
                p = [p p1];
                label = [label, {['TC manifold', num2str(i)]}];
            end
        end
        if (~isempty(xplot))&&(~isempty(xcolor))
            xplot = xplot(part,:);
            if ~input.Results.scatter_only
                dt = round(diff(tgrid(part)'));
                d = [0 find(dt>1) numel(part)];
                color = input.Results.line_color;
                for i=1:numel(d)-1
                    p1=plot(xplot(d(i)+1:d(i+1),1),xplot(d(i)+1:d(i+1),2),'Color',color);
                end
                p = [p p1];
                label = [label, {'trajectory'}];
            end
            if ~input.Results.line_only
                p1=scatter(xplot(:,1),xplot(:,2),50,xcolor(part),'Marker','.');
                p = [p p1];
            end
            label = [label, {'projected input'}];
        end
        if input.Results.legend
            legend(p,label);
        end
    case 3
        p = [];
        label = [];
        if ~isempty(xbackground)
            n_xbg = numel(xbackground);
            for i=1:n_xbg
                p1 = scatter3(xbackground{i}(:,1),xbackground{i}(:,2),xbackground{i}(:,3),50,'MarkerEdgeColor',[.8,.8,.8]*i/n_xbg,'Marker','.');
                label = [label, {['TC manifold', num2str(i)]}];
                p = [p p1];
            end
        end
        if (~isempty(xplot))&&(~isempty(xcolor))
            if ~input.Results.scatter_only
                dt = round(diff(tgrid(part)'));
                d = [0 find(dt>1) numel(part)]-1+part(1);
                color = input.Results.line_color;
                for i=1:numel(d)-1
                    p1=plot3(xplot(d(i)+1:d(i+1),1),xplot(d(i)+1:d(i+1),2),xplot(d(i)+1:d(i+1),3),'Color',color,'LineWidth',2);
                end
                p = [p,p1];
                label = [label,{'trajectory'}];
            end
            if ~input.Results.line_only
                p1=scatter3(xplot(part,1),xplot(part,2),xplot(part,3),120,xcolor(part),'Marker','.');
                p = [p,p1];
                label = [label,{'$\hat{x}$'}];
            end
            legend(p,label);
        end
        if input.Results.legend && ~isempty(xbackground)
            legend(p,label);
        end
        grid on
end
end

