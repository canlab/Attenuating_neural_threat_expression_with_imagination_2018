function ie_getstats = ie_getstats(IE, SE, NE, Plot_Title, Y_Label)
% This function gets group stats for imagined extinction
% Have the data for IE, SE, and NE in workspace

% quick check
if length(IE) ~= 20;warning('IE wrong size, check N');end;
if length(NE) ~= 24;warning('NE wrong size, check N');end;
if length(SE) ~= 22;warning('SE wrong size, check N');end;

% figure size
fig_pos_size=[.01 .3 .3 .3]; %top left corner, small square

% run plot and stats
pe_means=[mean(IE) mean(SE) mean(NE)];
pe_ste=[ste(IE) ste(SE) ste(NE)];

figure('units','normalized','position',fig_pos_size);
x=1:length(pe_means);wd = 0.6;
bar(x,[pe_means(1) nan nan],wd,'FaceColor',([80/255 162/255 206/255]),'EdgeColor',([0 0 0]),'LineWidth',1.5);
hold on
bar(x,[nan pe_means(2) nan],wd,'FaceColor',([179/255 179/255 179/255]),'EdgeColor',([0 0 0]),'LineWidth',1.5);
hold on
bar(x,[nan nan pe_means(3)],wd,'FaceColor',([76/255 76/255 76/255]),'EdgeColor',([0 0 0]),'LineWidth',1.5);
errorbar(x,pe_means,pe_ste,'.','Color','k','LineWidth',1.5);
ylabel(sprintf('%s',Y_Label),'FontSize',15);
% ylim([0 80]);
set(gca,'XTickLabel',{'IE','SE','NE'},'FontSize',13);
title(sprintf('%s',Plot_Title),'FontSize',14);

% run stats
% effect size: hedge's g here similar to cohens d
% most common measure of effect size for a One-Way ANOVA is Eta-squared (SSbtwn/SStot)
% Y: columns are groups, rows instances
fprintf('-------------------- output --------------------\n');
Y=NaN(24,3);
Y(1:size(IE,1),1)=IE;
Y(1:size(SE,1),2)=SE;
Y(1:size(NE,1),3)=NE;
[p,tbl,stats] = anova1(Y);
ES=cell2mat(tbl(2,2))/cell2mat(tbl(4,2));
F=cell2mat(tbl(2,5));
fprintf('---- %s----\n',Plot_Title);
if p < .05;disp('sig anova');end
fprintf('F(%d,%d)=%f, p=%f, Eta Sq=%f, means IE: %f SE: %f NE: %f\n',cell2mat(tbl(2,3)),cell2mat(tbl(3,3)),F,p,ES,stats.means(1),stats.means(2),stats.means(3));

% paired t-tests
[H,P,CI,STATS] = ttest2(Y(:,1),Y(:,2),'vartype','unequal');
ES=mes(Y(:,1),Y(:,2),'hedgesg');
if H > 1; disp('sig ttest');end
fprintf('IE v SE: t(%f)=%f,p=%f, CI=[%f %f], STD=[%f %f], Hedges g=%f\n',STATS.df,STATS.tstat,P,CI(1),CI(2),STATS.sd(1),STATS.sd(2),ES.hedgesg);
[H,P,CI,STATS] = ttest2(Y(:,2),Y(:,3),'vartype','unequal');
ES=mes(Y(:,2),Y(:,3),'hedgesg');
if H > 1; disp('sig ttest');end
fprintf('SE v NE: t(%f)=%f,p=%f, CI=[%f %f], STD=[%f %f], Hedges g=%f\n',STATS.df,STATS.tstat,P,CI(1),CI(2),STATS.sd(1),STATS.sd(2),ES.hedgesg);
[H,P,CI,STATS] = ttest2(Y(:,1),Y(:,3),'vartype','unequal');
ES=mes(Y(:,1),Y(:,3),'hedgesg');
if H > 1; disp('sig ttest');end
fprintf('IE v NE: t(%f)=%f,p=%f, CI=[%f %f], STD=[%f %f], Hedges g=%f\n',STATS.df,STATS.tstat,P,CI(1),CI(2),STATS.sd(1),STATS.sd(2),ES.hedgesg);


end