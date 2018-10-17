%% Attenuating neural threat expression with imagination
%  Reddan, Wager, & Schiller, Neuron 2018
% ------------------------------------------------------------
% this script uses code from
% https://github.com/canlab
% SPM8
% ------------------------------------------------------------
%% Supplementary Figure 2. Univariate activation during threat Acquisition 
% (CS + > CS-). This the FDR corrected contrast map (CS+ > CS-) resulting 
% from a univariate general linear model applied to all subjects (N=68) during 
% the Threat Acquisition phase. This map is similar to the multivariate Neural 
% Threat-Predictive Pattern in the main text. Regions with positive activations 
% include the dACC, PAG, amygdala and insula. Regions with negative activations 
% include the precuneus and orbitofrontal cortex. See Table S2 for a detailed 
% list of regions.
load IE_SuppFigure2_dat
cl=region(statsimg)
o2 = canlab_results_fmridisplay_marianne([],'full');
o2=addblobs(o2,cl,'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]});
o2=addblobs(o2,cl,'color',[0 0 0],'outline','transparent');
legend(o2);
cluster_table(cl,1,1,'writefile','IE_Supp_UnivariateAcqMap_Fdr05_k_cluster_table');

%% Supplementary 3 - Amygdala ROI
% load threat_predictice_weight_map
load IE_Neuron_SupplementaryFigure3_dat
%%group stats on CV PE
IE=PE_dat{1}(:,1)-PE_dat{1}(:,2);
SE=PE_dat{2}(:,1)-PE_dat{2}(:,2);
NE=PE_dat{3}(:,1)-PE_dat{3}(:,2);
Plot_Title='Amygdala Threat Pattern Expression during Recovery Test';
Y_Label='Average Pattern Expression';
ie_getstats(IE, SE, NE, Plot_Title, Y_Label);
ie_getstats_violin(IE, SE, NE, Plot_Title, Y_Label);

%% Supplementary 4 - Univariate late re-extinction
load IE_Neuron_SupplementaryFigure4_dat

% make 3 group whole brain thresholded plots of contrast
% one sample t-test, nothing survives fdr
IE_out =ttest(IE, .01, 'unc');
SE_out = ttest(SE, .01, 'unc');
NE_out = ttest(NE, .01, 'unc');

o2=canlab_results_fmridisplay_marianne([],'montagetype','full_long');
[o2, IEdat, IEsig, IEpcl, IEncl] = multi_threshold(IE_out, 'o2', o2, 'thresh', [.001 .005 .01], 'sizethresh', [5 5 5],'writestats');pause;removeblobs(o2);
[o2, SEdat, SEsig, SEpcl, SEncl] = multi_threshold(SE_out, 'o2', o2, 'thresh', [.001 .005 .01], 'sizethresh', [5 5 5],'writestats');pause;removeblobs(o2);
[o2, NEdat, NEsig, NEpcl, NEncl] = multi_threshold(NE_out, 'o2', o2, 'thresh', [.001 .005 .01], 'sizethresh', [5 5 5],'writestats');pause;removeblobs(o2);

cluster_table(region(fmri_data(which('IE_MultiThresh_T_P0_01_k5.nii'))), 1, 1,'writefile','LateReextCR_IE_Pos_Unc01_k5_cluster_table');
cluster_table(region(fmri_data(which('SE_MultiThresh_T_P0_01_k5.nii'))), 1, 1,'writefile','LateReextCR_SE_Pos_Unc01_k5_cluster_table');
cluster_table(region(fmri_data(which('NE_MultiThresh_T_P0_01_k5.nii'))), 1, 1,'writefile','LateReextCR_NE_Pos_Unc01_k5_cluster_table');
masks={'/rois/Hippocampus_CA1.nii';'/rois/bilateral_amygdala_50.ni'};
for r=1:length(masks)
    mask_name=strsplit(masks{r},'/');
    mask_name=mask_name{end}(1:end-4);
    
    dat=apply_mask(IE,masks{r});
    IEroi=mean(dat.dat)';
    AvgROI.(mask_name).IE=IEroi;
    
    dat=apply_mask(SE,masks{r});
    SEroi=mean(dat.dat)';
    AvgROI.(mask_name).SE=SEroi;
    
    dat=apply_mask(NE,masks{r});
    NEroi=mean(dat.dat)';
    AvgROI.(mask_name).NE=NEroi;
    
    Plot_Title=sprintf('%s',mask_name);
    Y_Label=sprintf('Threat Recovery\nAverage Univariate Beta Weight');
    ie_getstats(IEroi, SEroi, NEroi, Plot_Title, Y_Label);
    pause;
end

%% Supplementary 5 - SCR
% done in R - LMER
% IE_LMER_SCR_66.csv
% IE_LMER_SCR_42.csv
% IE_SuppFig4_MixedModel_SCR.R

%% Supplementary 5 - vmPFC pattern during extincton tp 3
% done in R
% IE_SupplementaryFigure5_dat_iesene_vmpfc_3rdtp_forpermute.csv
% IE_SuppFig5_Permute_Test_vmPFCspatialpattern_R.R