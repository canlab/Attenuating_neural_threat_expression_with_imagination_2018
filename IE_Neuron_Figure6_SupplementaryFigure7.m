%% Attenuating neural threat expression with imagination
%  Reddan, Wager, & Schiller, Neuron 2018
% ------------------------------------------------------------
% this script uses code from
% https://github.com/canlab
% SPM8
% ------------------------------------------------------------
% Figure 6. Brain regions that predict the success of imagined and real extinction. 
% A linear regression was performed on whole brain maps of differential activation 
% (CS+ > CS-) during extinction (time point 3) with participants? threat-predictive 
% pattern expression scores during the recovery test as the predictor. Because we 
% were most interested in brain regions that predicted extinction success, we 
% focused on the negatively correlated results (see Table S7 for a complete table 
% of activation). When thresholded and corrected for multiple comparisons (FDR corrected, 
% P < 0.05), both imagined (A) and real extinction (B) revealed activity in the vmPFC, 
% primary auditory cortex, and amygdala were related to extinction success. No significant
% correlates of extinction success were found in the no extinction group. C. Conjunction. 
% Signal in a priori ROIs that survived correction in the whole brain analysis was assessed. 
% The NAc predicted extinction success uniquely in the imagined extinction group, while 
% CA1 predicted extinction success uniquely in the real extinction group. Effect sizes 
% in the ROIs were estimated via max partial ?2 which ranges from 0 to 1. 

%% Run the OLS
load IE_Figure6_dat
gray_mask = fmri_data(which('gray_matter_mask.img'));
IE=apply_mask(fmri_data(filenames('/extinct_maps/IE*EXExt*.nii')),gray_mask);
SE=apply_mask(fmri_data(filenames('/extinct_maps/IE*NCExt*.nii')),gray_mask);
NE=apply_mask(fmri_data(filenames('/extinct_maps/IE*ACExt*.nii')),gray_mask);
% extract white matter
ie_m = extract_gray_white_csf(IE);
se_m = extract_gray_white_csf(SE);
ne_m = extract_gray_white_csf(NE);
% mean-center success scores and attach them to dat in dat.X
IE.X = [scale(IE_Y, 1) ie_m(:,2:3)];
SE.X = [scale(SE_Y, 1) se_m(:,2:3)];
NE.X = [scale(NE_Y, 1) ne_m(:,2:3)];
% % run the regresion, OLS with nuisnace
IE_out = regress(IE,.05, 'fdr');
SE_out = regress(SE,.05, 'fdr');
NE_out = regress(NE,.05, 'fdr');

%% A - Imagined Extinction
IEdat=IE_out.t;ie_fdr_thresh=IE_out.t.threshold(1); % T
IEdat.thr_type='FDR';IEdat.threshold=ie_fdr_thresh;
IEdat.p(:,2:end)=[];
IEdat.ste(:,2:end)=[];
IEdat.sig(:,2:end)=[];
IEdat.dat(:,2:end)=[];
IEdat.dat=IEdat.dat .* (IEdat.dat<0);
IEdat.p=IEdat.p .* (IEdat.dat<0);
IEdat.ste=IEdat.ste .* (IEdat.dat<0);
IEdat.sig=IEdat.sig .* (IEdat.dat<0);
IEdat.dat=IEdat.dat .* (IEdat.sig>0);
cluster_table(region(IEdat),0,1,'writefile','IE_3rdExtCorrThreatPE_wmcontrol_OLS_fdr05_k1_tstat');
o2=canlab_results_fmridisplay_marianne([],'montagetype','full');
o2= addblobs(o2,region(IEdat));
saveas(gcf, 'IE_3rdExtCorrThreatPE_wmcontrol_OLS_fdr05_k1_tmap.png');
legend(o2);saveas(gcf, 'IE_3rdExtCorrThreatPE_wmcontrol_OLS_fdr05_k1_tmap.png');
close;
%% B - Real Extinction 
SEdat=SE_out.t;se_fdr_thresh=SE_out.t.threshold(1); % T
SEdat.thr_type='FDR';SEdat.threshold=se_fdr_thresh;
SEdat.p(:,2:end)=[];
SEdat.ste(:,2:end)=[];
SEdat.sig(:,2:end)=[];
SEdat.dat(:,2:end)=[];
SEdat.dat=SEdat.dat .* (SEdat.dat<0);
SEdat.p=SEdat.p .* (SEdat.dat<0);
SEdat.ste=SEdat.ste .* (SEdat.dat<0);
SEdat.sig=SEdat.sig .* (SEdat.dat<0);
SEdat.dat=SEdat.dat .* (SEdat.sig>0);
cluster_table(region(SEdat),0,1,'writefile','SE_3rdExtCorrThreatPE_wmcontrol_OLS_fdr05_k1_b_tstat');
o2=canlab_results_fmridisplay_marianne([],'montagetype','full');
o2= addblobs(o2,region(SEdat));
saveas(gcf, 'SE_3rdExtCorrThreatPE_wmcontrol_OLS_fdr05_k1_tmap.png');
legend(o2);saveas(gcf, 'SE_3rdExtCorrThreatPE_wmcontrol_OLS_fdr05_k1_tmap.png');
close;

%% Alt ways of getting montages
title='SuccessRegress_Neg_OLS_map_k5_fdr_';

% slices to get
% coronal, y = = -5 (amyg) & +10 (NAc)
% sag, x = -3 (vmpfc), + 21 (CA1)
out_col=[.75 .75 .75];
for g=1:2
    if g==1
        cl=region(IEdat);
        name='IE_neg_overlap';
    elseif g==2
        cl=region(SEdat);
        name='SE_neg_overlap';
    end

o2 = fmridisplay;o2 = montage(o2, 'saggital', 'slice_range', [-4 -2], 'onerow', 'spacing', 1);
o2= addblobs_kthresh(o2,cl, 'sizethresh', 5);
o2=addblobs(o2,cl,'outline','outline_color',out_col);
saveas(gcf, sprintf('%s_sag_n4ton2_%s.png',title, name));
% 
close all;
o2 = fmridisplay;o2 = montage(o2, 'saggital', 'slice_range', [21 23], 'onerow', 'spacing', 1);
o2= addblobs_kthresh(o2,cl, 'sizethresh', 5);
o2=addblobs(o2,cl,'outline','outline_color',out_col);
saveas(gcf, sprintf('%s_sag_21top23_%s.png',title, name));

close all;
o2 = fmridisplay;o2 = montage(o2, 'coronal', 'slice_range', [-5 -3], 'onerow', 'spacing', 1);
o2= addblobs_kthresh(o2,cl, 'sizethresh', 5);
o2=addblobs(o2,cl,'outline','outline_color',out_col);
saveas(gcf, sprintf('%s_coron_n5ton3_%s.png',title, name));

close all;
o2 = fmridisplay;o2 = montage(o2, 'coronal', 'slice_range', [9 11], 'onerow', 'spacing', 1);
o2= addblobs_kthresh(o2,cl, 'sizethresh', 5);
o2=addblobs(o2,cl,'outline','outline_color',out_col);
saveas(gcf, sprintf('%s_coron_9to11_%s.png',title, name));

legend(o2);
saveas(gcf, sprintf('%s_LEGEND_%s.png',title, name));

close all;
end

%% C - Conjunction

% get unthresholded neg vals
IEdat_nothresh=IE_out.t;
IEdat_nothresh.p=[];
IEdat_nothresh.ste(:,2:end)=[];
IEdat_nothresh.sig=[];
IEdat_nothresh.dat(:,2:end)=[];
IEdat_nothresh.dat=IEdat_nothresh.dat .* (IEdat_nothresh.dat<0);

SEdat_nothresh=SE_out.t;
SEdat_nothresh.p=[];
SEdat_nothresh.ste(:,2:end)=[];
SEdat_nothresh.sig=[];
SEdat_nothresh.dat(:,2:end)=[];
SEdat_nothresh.dat=SEdat_nothresh.dat .* (SEdat_nothresh.dat<0);


% load mask objects
amyg=fmri_data('/Users/maus/Desktop/2018_NN_Paper_Dat_and_Stat/IE_Figure4D/rois/bilateral_amygdala_50.nii');
vm=fmri_data('/Users/maus/Desktop/2018_NN_Paper_Dat_and_Stat/IE_Figure4D/rois/vmpfc_50perc_binary.nii');
aud=fmri_data('/Users/maus/Desktop/2018_NN_Paper_Dat_and_Stat/IE_Figure4D/rois/Auditory_Te10.nii');
nac=fmri_data('/Users/maus/Desktop/2018_NN_Paper_Dat_and_Stat/IE_Figure4D/rois/Accumbens_mask.nii');
ca1=fmri_data('/Users/maus/Desktop/2018_NN_Paper_Dat_and_Stat/IE_Figure4D/rois/Hippocampus_CA1.nii');

masks=filenames('/Users/maus/Desktop/2018_NN_Paper_Dat_and_Stat/IE_Figure6/rois/*.nii', 'absolute');
% make the contour maps of the ROIs
% apply to unthresh, IE_out.b // or t

view=[2,3,1,2,1]; % which view for mask#
coord=[10,6,21,-10,-2]; % which coord for mask#
for m=1:length(masks)
    mask_name=strsplit(masks{m},'/');
    mask_name=mask_name{end}(1:end-4);
    % unthresh    
    se_roi=apply_mask(fmri_data(SEdat_nothresh_all), fmri_data(masks{m}));
    ie_roi=apply_mask(fmri_data(IEdat_nothresh_all), fmri_data(masks{m}));
    
    info = roi_contour_map([region(ie_roi) region(se_roi)], 'cluster', 'use_same_range', 'colorbar','xyz',view(m),'xyz_coord',coord(m));
    
    masked_dat.mask_name{m}=mask_name;
    masked_dat.info{m}=info;
    masked_dat.min{m}=[min(abs(info{1}.Z(find(info{1}.Z<0)))) min(abs(info{2}.Z(find(info{2}.Z<0)))) min(abs(info{3}.Z(find(info{3}.Z<0)))) min(abs(info{4}.Z(find(info{4}.Z<0))))];
    masked_dat.max{m}=[min((info{1}.Z)) min((info{2}.Z)) min((info{3}.Z)) min((info{4}.Z))];
    masked_dat.notes='first col IE, second col SE';
    
    masked_dat.min{m}
    masked_dat.max{m}
    pause;
     saveas(gcf, sprintf('SucReg_unthresh_contourplot_t_%s.png',mask_name)); 
    
    clear ie_roi;
    clear se_roi;
    clear info;
 
end

%% Get effect sizes for each voxel
IE=apply_mask(fmri_data(filenames('/extinct_maps/IE*EXExt*.nii')),gray_mask);
SE=apply_mask(fmri_data(filenames('/extinct_maps/IE*NCExt*.nii')),gray_mask);
% extract white matter
ie_m = extract_gray_white_csf(IE);
se_m = extract_gray_white_csf(SE);
% mean-center success scores and attach them to dat in dat.X
IE.X = [scale(IE_Y, 1) ie_m(:,2:3)];
SE.X = [scale(SE_Y, 1) se_m(:,2:3)];

for g=1:2
    if g==1
        dat=IE;IEes=IEdat;IEes.fullpath=[];IEes.dat_descrip='Partial eta sq for Threat Pred';IEes.dat=[];
    elseif g ==2
        dat=SE;SEes=SEdat;SEes.fullpath=[];SEes.dat_descrip='Partial eta sq for Threat Pred';SEes.dat=[];
    end
    
    X = intercept(dat.X, 'add');
    for i = 1:size(dat.dat,1)
        [B,~,~,~,~] = regress(dat.dat(i,:)',X);
        
        y_mean=mean(dat.dat(i,:)');
        y_hat=y_mean+(B(1)*X(:,1)); % just care about first predictor, threat
        SST=sumsqr(dat.dat(i,:)'-y_mean); %total sum of squares
        SSR=sumsqr(y_hat - y_mean); %explained sum of squares
        SSE=sumsqr(y_hat - dat.dat(i,:)'); %residual Y actual from Y predicted
        ParEtaSq(i)=SSR/(SSE+SSR);
        
        %     % don't need intercept
        %     tbl = table(X(:,1),X(:,2),X(:,3),dat.dat(i,:)','VariableNames',{'Threat','WhiteMatter','CSF','Brain'});
        %     lm = fitlm(tbl,'Brain~Threat+WhiteMatter+CSF')
        %     PartialEtaSq=lm.SSR/(lm.SSR+lm.SSE); % ssr = for the effect, sse = err
        
    end
    if g==1
        IEes.dat=ParEtaSq';
    elseif g==2
        SEes.dat=ParEtaSq';
    end
    clear ParEtaSq;
end

for m=1:length(masks)
    mask_name=strsplit(masks{m},'/');
    mask_name=mask_name{end}(1:end-4)
    % thresh    
    seroi=apply_mask(fmri_data(SEes),fmri_data(masks{m}));
%     plot(seroi);
    secl=region(seroi);
%     for x=1:length(secl)
%     fprintf('SE %s: NumVox=%d, avgES=%f, steES=%f, maxES=%f, coor=[%d %d %d]',mask_name, size(seroi.dat,1),mean(seroi.dat,1),ste(seroi.dat),max(seroi.dat),secl(x).mm_center(1),secl(x).mm_center(2),secl(x).mm_center(3));
%     pause;
%     end
    ROI.max_es(m,2)=max(seroi.dat);
    ROI.mean_es(m,2)=mean(seroi.dat,1);
    ROI.ste_es(m,2)=ste(seroi.dat);
    
    ieroi=apply_mask(fmri_data(IEes),fmri_data(masks{m}));
%     plot(ieroi);
    iecl=region(ieroi);
%     for x=1:length(iecl)
%         fprintf('IE %s: NumVox=%d, avgES=%f, steES=%f, maxES=%f, coor=[%d %d %d]',mask_name, size(ieroi.dat,1),mean(ieroi.dat,1),ste(ieroi.dat),max(ieroi.dat),iecl(x).mm_center(1),iecl(x).mm_center(2),iecl(x).mm_center(3));
%     pause;
%     end
    ROI.max_es(m,1)=max(ieroi.dat);
    ROI.mean_es(m,1)=mean(ieroi.dat,1);
    ROI.ste_es(m,1)=ste(ieroi.dat);
    
end

%% Supplementary Figure 7
% Supplementary Figure 7. Distribution of data supporting the effects in ROIs 
% that predict the success of imagined and real extinction. Scatterplots are 
% shown for descriptive purposes, and illustrate the distribution of individual 
% data values in Figure 6, but should not be taken as indicative of the true 
% effect sizes [71]. 

ie_color_blue=[80/255 162/255 206/255];
se_color_gray=[179/255 179/255 179/255];
for m=3:length(masks)
    mask_name=strsplit(masks{m},'/');
    mask_name=mask_name{end}(1:end-4);
    
    % make biased mask of cluster
    ie_clust=apply_mask(fmri_data(IEdat),masks{m});ie_clust.dat=logical(ie_clust.dat);
    se_clust=apply_mask(fmri_data(SEdat),masks{m});se_clust.dat=logical(se_clust.dat);
    
    IE=apply_mask(fmri_data(filenames('/extinct_maps/IE*EXExt*.nii')),ie_clust);
    SE=apply_mask(fmri_data(filenames('/extinct_maps/IE*NCExt*.nii')),se_clust);
    ieY=mean(IE.dat)';
    seY=mean(SE.dat)';
    % mean-center success scores and add nuisance and intercept
    ieX = [scale(IE_Y, 1) ie_m(:,2:3) ones(length(IE_Y),1)];
    seX = [scale(SE_Y, 1) se_m(:,2:3) ones(length(SE_Y),1)];
    
    [ieB,ieBINT,ieR,ieRINT,ieSTATS] = regress(ieY,ieX);iep=ieSTATS(3);
    [seB,seBINT,seR,seRINT,seSTATS] = regress(seY,seX);sep=seSTATS(3);
    
    ieParEtaSq=get_ParEtaSq(ieY,ieB(1),ieX(:,1));
    seParEtaSq=get_ParEtaSq(seY,seB(1),seX(:,1));
    
    set(gcf,'units','normalized','position',[.1 .5 .2 .2]);
    scatter(ieY,ieX(:,1),30,ie_color_blue,'filled'); grid on; hold on;
    plot(ieX(:,1),ieB(4)+ieB(1)*ieX(:,1),'b','LineWidth',1);
    
    scatter(seY,seX(:,1),30,se_color_gray,'filled'); grid on; hold on;
    plot(seX(:,1),seB(4)+seB(1)*seX(:,1),'k','LineWidth',1);
    
    DblDip{m}.mask_name=mask_name;
    DblDip{m}.b=[ieB,seB];
    DblDip{m}.p=[iep,sep];
    DblDip{m}.paretasq=[ieParEtaSq,seParEtaSq];
    DblDip{m}.rsquare=[ieSTATS(1),seSTATS(1)];
    DblDip{m}.notes='first col IE, second SE, unbiased ROI avg bet fig 6, plot name Fig6_DoubleDipPlot_';
    title(sprintf('%s',mask_name),'FontSize',16);
%     xlabel('Threat Expression during Recovery Test','FontSize',16);
%     ylabel(sprintf('Average Differential Activation (CS+ > CS-) \n during Extinction Timepoint 3'),'FontSize',16);
    ylabel(sprintf('Average Beta'),'FontSize',16);
    xlabel('Threat Recovery','FontSize',16);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    saveas(gcf, sprintf('Fig6_DoubleDipPlot_%s.png',mask_name));
    display(sprintf('in the bilateral %s of the Imagined Extinction group (b = %.2f,p = %.2f,partial = %.2f Rsq = %.2f', mask_name,ieB(1),iep(1),ieParEtaSq,ieSTATS(1)));
    display(sprintf('in the bilateral %s of the Real Extinction group (b = %.2f,p = %.2f,partial = %.2f Rsq = %.2f', mask_name,seB(1),sep(1),seParEtaSq,seSTATS(1)));

    pause;close all;
%     format short;
end
save Fig6_CorrPlots_BiasedCluster DblDip