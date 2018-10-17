function ie_get_lineplot = ie_get_lineplot(IE, SE, NE, Plot_Title, Y_Label)
% This function takes the data to make the extinction time course plot
% Subj-wise observations in rows, Timepoints in Cols

% ste is from canlab CanlabCore/Statistics_tools/ste.m
% nanmean is matlab built in not field trip

% quick check
if length(IE) ~= 20;warning('IE wrong size, check N');end;
if length(NE) ~= 24;warning('NE wrong size, check N');end;
if length(SE) ~= 22;warning('SE wrong size, check N');end;

% colors ... 
% IE 50A2CE / [80/255 162/255 206/255]
% SE B3B3B3 / [179/255 179/255 179/255]
% NE 4C4C4C / [76/255 76/255 76/255]
fs=18;
sz=8;
lnsz=2;

for g=1:3
    if g==1;
        name='Imagined Extinction';%pos=0;
        Y=mean(IE);err=ste(IE);
        cf=[80/255 162/255 206/255];
    elseif g==2
        name='Real Extinction';%pos=3;
        Y=mean(SE);err=ste(SE);
        cf=[179/255 179/255 179/255];
    elseif g==3
        name='No Extinction';%pos=6;
        Y=mean(NE);err=ste(NE);
        cf=[76/255 76/255 76/255];
    end
    x=[1:1:size(Y,2)];
    errorbar(x,Y,err,'-o','LineWidth',lnsz,'MarkerSize',sz,...
    'MarkerEdgeColor',cf-.1,'MarkerFaceColor',cf,'Color',cf);hold on;
end
ylabel(Y_Label,'FontSize',fs);
title(sprintf('%s',Plot_Title),'FontSize',fs+2);
xt = get(gca, 'XTick');
set(gca, 'FontSize', fs)
r=refline(0,0);
r.Color='k';
r.LineWidth=.35;
legend('Imagined Extinction','Real Extinction','None', 'Location', 'southeast');
saveas(gcf,sprintf('%s_%s_TimeCoursePlot',Plot_Title,Y_Label),'fig');
saveas(gcf,sprintf('%s%s__TimeCoursePlot',Plot_Title,Y_Label),'png');

end

