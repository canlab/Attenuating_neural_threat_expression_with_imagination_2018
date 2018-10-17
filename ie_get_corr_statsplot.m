function ie_get_corr_statsplot = ie_get_corr_statsplot(IE, SE, NE, Plot_Title, Y_Label, X_Label)
% This function gets correlation stats for imagined extinction study
% Have the data for IE, SE, and NE in workspace so that col 1 is X and col
% 2 is Y

% quick check
if length(IE) ~= 20;warning('IE wrong size, check N');end;
if length(NE) ~= 24;warning('NE wrong size, check N');end;
if length(SE) ~= 22;warning('SE wrong size, check N');end;

% colors ... 
% IE 50A2CE / [80/255 162/255 206/255]
% SE B3B3B3 / [179/255 179/255 179/255]
% NE 4C4C4C / [76/255 76/255 76/255]
fs=18;
sz=70;
lnsz=1.5;
ie_cf=[80/255 162/255 206/255];
se_cf=[179/255 179/255 179/255];
ne_cf=[76/255 76/255 76/255];
% all_cf=[repmat(ie_cf,20,1);repmat(se_cf,22,1);repmat(ne_cf,24,1)];

for g=1:3
    if g==1;
        name='Imagined Extinction';cf=ie_cf;%pos=0;
        X=IE(:,1);Y=IE(:,2);
    elseif g==2
        name='Real Extinction';cf=se_cf;%pos=3;
        X=SE(:,1);Y=SE(:,2);
    elseif g==3
        name='No Extinction';cf=ne_cf;%pos=6;
        X=NE(:,1);Y=NE(:,2);
    end
        set(gcf,'units','normalized','position',[.1 .5 .9 .5]);
%         pos=pos+1;
        subplot(1,3,g);
        s=scatter(X,Y,sz,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',cf,'LineWidth',lnsz);
        [r,p]=corr(X,Y);
        % do OLS        
        s=lsline;
        if p > 0.05
            set(s,'LineStyle','--','color','k');
        else
            set(s,'LineStyle','-','color','k');
        end
        [r,p] = corr(X,Y);
        disp(sprintf('%s Pearsons Correlation, %s, R=%f, p=%f', name, Plot_Title, r, p));
        
%         % do robust
%         [brob STATS] = robustfit(X,Y);
%         hold on;
%         if p > 0.05
%             plot(X,brob(1)+brob(2)*X,'--k','LineWidth',lnsz);
%         else
%             plot(X,brob(1)+brob(2)*X,'-k','LineWidth',lnsz)
%         end

        
        xlabel(X_Label,'FontSize',fs);
        ylabel(Y_Label,'FontSize',fs);
        title(sprintf('%s \n %s \n R=%g,p=%g',Plot_Title,name,r,p),'FontSize',fs+2);
        linkaxes;
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', fs);
end
saveas(gcf,sprintf('%s_Corr_%s_%s',Plot_Title,X_Label,Y_Label),'fig');
saveas(gcf,sprintf('%s_Corr_%s_%s',Plot_Title,X_Label,Y_Label),'png');