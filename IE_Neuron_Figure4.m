%% Attenuating neural threat expression with imagination
%  Reddan, Wager, & Schiller, Neuron 2018
% ------------------------------------------------------------
% this script uses code from
% https://github.com/canlab
% SPM8
% ------------------------------------------------------------
% Figure 4. Neural Networks Supporting Extinction Learning. Functional 
% connectivity between a priori regions of interest (vmPFC, hippocampus 
% (CA1), laterobasal amygdala (Amg-LB), central nucleus of the amygdala 
% (Amyg-CM), NAc, PAG, and the primary auditory cortex) across the entire 
% extinction session was investigated without a contrast. Links are plotted 
% if the partial correlation test between two nodes was significant after 
% FDR correction within each group. The link length between each pair 
% of nodes is related to the absolute value of the t-statistic associated 
% with the partial correlation between each pair, however, link length is a 
% graphical representation of this value, adjusted to exist in a 2D space. Here, 
% shorter links represent greater functional connectivity. The size of each node 
% represents betweenness centrality, which is the number of shortest paths that 
% pass through the node, rescaled for plotting purposes with a sigmoid function. 
% Node colors indicate clusters determined by ward linkage. A. Imagined Extinction 
% Connectivity Network. The most central nodes (i.e., nodes most connected with 
% other nodes) of the imagined extinction network included the right Nac, vmPFC, 
% left Amg-CM, and right Primary Auditory Cortex. B. Real Extinction Connectivity 
% Network. The most central nodes of real extinction included the left CA1, vmPFC, 
% and right Nac. C. No Extinction Connectivity Network. The most central nodes 
% during the no extinction, unrelated imagination included left NAc, left Amyg-CM, 
% and right NAc. D. Comparison of Nodal Centrality Across Networks. 10 out of 12 
% nodes yielded group significant differences in betweenness centrality in a one-way 
% ANOVA test, FDR corrected for multiple comparisons (p < 0.05). Asterisks represent 
% pairwise significant difference via a post hoc t-test. The most central nodes to 
% imagined extinction are indicated by a dashed box. E. Each Network is Unique.  In 
% order to test the separability of each matrix, and thereby assess their relative 
% similarity, we trained three different binary linear support vector machine classifiers 
% on participant connectivity data (the partial correlations between each pair of nodes). 
% Leave-one-out cross-validated accuracy scores revealed that each group was linearly 
% separable from one another with high accuracy (Real v None:  91% +- 4.2% STE; 
% Imagined v None: 89% +- 4.8% STE; Imagined v Real: 88% +- 5.0% STE). 

%% Connectivity Analysis
% IE_MultiVar_Connect_Final_NoContrast_2017.m

%% A-C force directed graphs / Connectivity Networks.

load IE_Figure4_dat
% addpath('/CANLabRepos/trunk/densityUtility/additional_utilities_and_imgs/')

% OUT came from: canlab_connectivity_predict(dat, subject_grouping, 'partialr');
% t's from: OUT.stats.fdr_thresholded_tvalues;

% change blending options of patch faces in the editor to get rid of the
% grey bands on the brain image

% A - Imagined Extinction Connectivity Network
OUT = ieOUT;
t = iet;
[ie_stats, handles] = canlab_force_directed_graph(t,'rset', OUT.parcelindx,'cl',cl);
close all;

% B - Real Extinction Connectivity Network
OUT = seOUT;
t = set;
[se_stats, handles] = canlab_force_directed_graph(t,'rset', OUT.parcelindx,'cl',cl);

% C - No Extinction Connectivity Network
OUT = neOUT;
t = net;
[ne_stats, handles] = canlab_force_directed_graph(t,'rset', OUT.parcelindx,'cl',cl);

%% D - Comparison of Node Centrality
for i = 1:size(ieOUT.betweenness_centrality,2)
    fprintf('Node: %s\n',names{i})
    % set up dat
    IE=ieOUT.betweenness_centrality(:,i);
    SE=seOUT.betweenness_centrality(:,i);
    NE=neOUT.betweenness_centrality(:,i);
    Plot_Title=names{i};
    Y_Label='Betweenness Centrality';
    ie_getstats(IE, SE, NE, Plot_Title, Y_Label);
    pause;
    figure(1);
    saveas(gcf,sprintf('nodesize_%s.png',Plot_Title));
    close all; 
end

% subplot
fig_pos_size=[.01 .3 .3 .9];figure('units','normalized','position',fig_pos_size);

for i = 1:size(ieOUT.betweenness_centrality,2)
    subplot(3,4,i)
    IE=ieOUT.betweenness_centrality(:,i);
    SE=seOUT.betweenness_centrality(:,i);
    NE=neOUT.betweenness_centrality(:,i);
    pe_means=[mean(IE) mean(SE) mean(NE)];
    pe_ste=[ste(IE) ste(SE) ste(NE)];
    x=1:length(pe_means);wd = 0.6;
    bar(x,[pe_means(1) nan nan],wd,'FaceColor',([80/255 162/255 206/255]),'EdgeColor',([0 0 0]),'LineWidth',1.5);
    hold on
    bar(x,[nan pe_means(2) nan],wd,'FaceColor',([179/255 179/255 179/255]),'EdgeColor',([0 0 0]),'LineWidth',1.5);
    hold on
    bar(x,[nan nan pe_means(3)],wd,'FaceColor',([76/255 76/255 76/255]),'EdgeColor',([0 0 0]),'LineWidth',1.5);
    errorbar(x,pe_means,pe_ste,'.','Color','k','LineWidth',1.5);
    ylim([0 80]);
    title(sprintf('%s',names{i}),'FontSize',14); 
    hold off;
end
saveas(gcf,sprintf('nodesize_subplot_allnodes.png'));

%% E - Each Network is Unique
load IE_Figure4E_Compare_IE_v_NE
ROC = roc_plot(OUT.PREDICT.pairwise_association.other_output{2}, logical(OUT.PREDICT.pairwise_association.Y > 0));
hold on;

load IE_Figure4E_Compare_IE_v_SE
ROC = roc_plot(OUT.PREDICT.pairwise_association.other_output{2}, logical(OUT.PREDICT.pairwise_association.Y > 0));
hold on;

load IE_Figure4E_Compare_SE_v_NE
ROC = roc_plot(OUT.PREDICT.pairwise_association.other_output{2}, logical(OUT.PREDICT.pairwise_association.Y > 0));

%  
% For pairwise
% ie v ne
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	-0.22	Sens:	 85% CI(71%-96%)	Spec:	 92% CI(81%-100%)	PPV:	 89% CI(75%-100%)	Nonparametric AUC:	0.95	Parametric d_a:	2.18	  Accuracy:	 89% +- 4.8% (SE), P = 0.000000
%  
% Ie v se
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	-0.16	Sens:	 85% CI(71%-100%)	Spec:	 91% CI(80%-100%)	PPV:	 89% CI(76%-100%)	Nonparametric AUC:	0.92	Parametric d_a:	1.90	  Accuracy:	 88% +- 5.0% (SE), P = 0.000000
%  
% Se v ne
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	-0.10	Sens:	100% CI(100%-100%)	Spec:	 83% CI(71%-96%)	PPV:	 85% CI(73%-96%)	Nonparametric AUC:	0.97	Parametric d_a:	2.50	  Accuracy:	 91% +- 4.2% (SE), P = 0.000000
