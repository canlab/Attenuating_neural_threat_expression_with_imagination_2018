%% Attenuating neural threat expression with imagination
%  Reddan, Wager, & Schiller, Neuron 2018
% ------------------------------------------------------------
% this script uses code from
% https://github.com/canlab
% SPM8
% ------------------------------------------------------------
% Figure 2. Multivariate Neural Threat-Predictive Pattern Trained on Acquisition 
% A. Neural Threat-Predictive Pattern. A distributed pattern of threat was developed 
% using a linear support vector machine trained on subject-wise (N=68) univariate 
% brain activation to the unreinforced CS+ (threatening) and CS- (non-threatening) 
% stimuli during threat acquisition. When thresholded (bootstrapped 5,000 samples) 
% and corrected for multiple comparisons (FDR p < 0.05), the ?threat-predictive 
% pattern? revealed a distributed network of threat representation including the 
% PFC, periaqueductal gray (PAG), insula, subregions of the basal ganglia, and dorsal 
% anterior cingulate (dACC). B. ROC Plot. The neural threat-predictive pattern 
% yielded a classification accuracy of 77% (+- 3.6 STE, p < .0001, AUC: 0.82) in 
% a modified leave-three-subjects-out cross-validation (CV) procedure. One subject 
% from each of the three groups was left out of each fold to reduce potential bias 
% from group assignment.

%% A - Neural Threat-Predicitve Pattern
load IE_Figure2A_dat
obj=statistic_image();
obj.dat=stats.WTS.wZ';
obj.p=stats.WTS.wP';
obj.volInfo=stats.weight_obj.volInfo;
obj = threshold(obj, .05, 'fdr', 'k', 1, 'mask', which('gray_matter_mask.img'))
cl=region(obj);
o2 = montage(o2, 'sagittal', 'slice_range', [-3 3], 'onerow', 'spacing', 1);
o2 = montage(o2, 'axial', 'slice_range', [-20 10], 'onerow', 'spacing', 5);
o2 = montage(o2, 'coronal', 'slice_range', [-6 0], 'onerow','spacing', 2);
o2=addblobs(o2,cl,'splitcolor', {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]})
o2=addblobs(o2,cl,'color',[0 0 0],'outline','transparent')
table=cluster_table(cl, 0, 0,'writefile','ThreatSig_FDR05_K0_clusttable');
select=[15 16 24 49];
info = roi_contour_map(cl(select), 'cluster', 'use_same_range', 'colorbar');

% obj.dat=obj.dat .* obj.sig;
obj.fullpath=['Figure2A_NeuralThreatPredictivePattern.nii'];
write(obj,'thresh','keepdt');


%% B - ROC Plot
load IE_Figure2B_dat
ROC = roc_plot(stats.dist_from_hyperplane_xval, logical(stats.Y > 0), 'unpaired');
% -------------------- output --------------------
% ROC_PLOT Output: Single-interval, Optimal overall accuracy
% Threshold:	0.21	Sens:	 71% CI(61%-79%)	Spec:	 84% CI(76%-91%)	PPV:	 81% CI(72%-90%)	Nonparametric AUC:	0.82	Parametric d_a:	1.22	  Accuracy:	 77% +- 3.6% (SE), P = 0.000000
