%% IMAGINED EXTINCTION - Multivariate Functional Connectivity Analysis
%   by marianne, 2017
% functions:
%   canlab_connectivity_predict
%   canlab_connectivity_preproc

%% STEP 1 - preprocessing the signal
% rm subjects: 270 (AC) 231 (NC)
cd('/Volumes/engram/Users/marianne/Documents/MATLAB/MultiVarConnect/');
sub_dirs=filenames('/Volumes/engram/labdata/data/Imagination/Imaging/IE*EX/','absolute');
subjects=filenames('/Volumes/engram/labdata/data/Imagination/Imaging/IE*EX/Functional/Preprocessed/r2EXT/swraIE*.nii','absolute');
numsubs=length(subjects);
covs=filenames('/Volumes/engram/labdata/data/Imagination/Imaging/IE*EX/Functional/Preprocessed/r2EXT/spm_modeling/noise_model_1.mat','absolute');
addpath(genpath('/Volumes/engram/Resources/Respository/trunk/densityUtility/'));

mask_dir='/Volumes/engram/Users/marianne/Documents/MATLAB/MultiVarConnect/ROIs/*.nii';
roi_masks=filenames([mask_dir, '*']);num_masks=size(roi_masks,1);
roi_masks(12)=[];
roi_masks(5)=[];
for r = 1:length(roi_masks)
    mask_split = strsplit(roi_masks{r},'/');
    mask_name = mask_split{end}(1:end-4);
    cl(r) = region(fmri_data(roi_masks{r}));
    names{r} = mask_name;
end

save IE_MultiVar_Connect_SetUp
%%
TR=2;
for n = 1:numsubs
    dat = fmri_data(subjects{n});
    subject_dir = sub_dirs{n};
    load(covs{n});
    dat.covariates=R;
    [preprocessed_dat, roi_val] = canlab_connectivity_preproc(dat, 'vw', 'linear_trend','datdir',subject_dir, 'bpf', [.008 .25], TR, 'extract_roi', roi_masks, 'average_over'); % OR 'unique_mask_values'
    ie_dat{n}.preprocessed=preprocessed_dat;
    ie_dat{n}.rois=roi_val;
    ie_dat{n}.subj=subjects{n};
end

save IE_MultiVarConnectivity_2017_Final_IE -v7.3 ie_dat

% rm subjects: 270 (AC) 231 (NC)
sub_dirs=filenames('/Volumes/engram/labdata/data/Imagination/Imaging/IE*NC/','absolute');
sub_dirs(8)=[];
subjects=filenames('/Volumes/engram/labdata/data/Imagination/Imaging/IE*NC/Functional/Preprocessed/r2EXT/swraIE*.nii','absolute');
subjects(8)=[];
numsubs=length(subjects);
covs=filenames('/Volumes/engram/labdata/data/Imagination/Imaging/IE*NC/Functional/Preprocessed/r2EXT/spm_modeling/noise_model_1.mat','absolute');
covs(8)=[];

for n = 1:numsubs
    dat = fmri_data(subjects{n});
    subject_dir = sub_dirs{n};
    load(covs{n});
    dat.covariates=R;
    [preprocessed_dat, roi_val] = canlab_connectivity_preproc(dat, 'vw', 'linear_trend','datdir',subject_dir, 'bpf', [.008 .25], TR, 'extract_roi', roi_masks, 'average_over');
%     [preprocessed_dat, roi_val] = canlab_connectivity_preproc(dat, 'vw', 'linear_trend','datdir',subject_dir, 'bpf', [.008 .25], TR, 'extract_roi', parcel_mask, 'unique_mask_values');
    se_dat{n}.preprocessed=preprocessed_dat;
    se_dat{n}.rois=roi_val;
    se_dat{n}.subj=subjects{n};
end

save IE_MultiVarConnectivity_2017_Final_SE -v7.3 se_dat

% rm subjects: 270 (AC) 231 (NC)
sub_dirs=filenames('/Volumes/engram/labdata/data/Imagination/Imaging/IE*AC/','absolute');
sub_dirs(13)=[];
subjects=filenames('/Volumes/engram/labdata/data/Imagination/Imaging/IE*AC/Functional/Preprocessed/r2EXT/swraIE*.nii','absolute');
subjects(13)=[];
numsubs=length(subjects);
covs=filenames('/Volumes/engram/labdata/data/Imagination/Imaging/IE*AC/Functional/Preprocessed/r2EXT/spm_modeling/noise_model_1.mat','absolute');
covs(13)=[];

for n = 1:numsubs
    dat = fmri_data(subjects{n});
    subject_dir = sub_dirs{n};
    load(covs{n});
    dat.covariates=R;
    [preprocessed_dat, roi_val] = canlab_connectivity_preproc(dat, 'vw', 'linear_trend','datdir',subject_dir, 'bpf', [.008 .25], TR, 'extract_roi', roi_masks, 'average_over');
%     [preprocessed_dat, roi_val] = canlab_connectivity_preproc(dat, 'vw', 'linear_trend','datdir',subject_dir, 'bpf', [.008 .25], TR, 'extract_roi', parcel_mask, 'unique_mask_values');
    ne_dat{n}.preprocessed=preprocessed_dat;
    ne_dat{n}.rois=roi_val;
    ne_dat{n}.subj=subjects{n};
end

save IE_MultiVarConnectivity_2017_NE -v7.3 ne_dat

%%
load IE_MultiVarConnectivity_2017_IE
load IE_MultiVarConnectivity_2017_SE
load IE_MultiVarConnectivity_2017_NE

%%
close all;
clear sub_rois;
a=1;
for s = 1:size(ie_dat,2)
   for r = 1:size(ie_dat{s}.rois,1)
        sub_rois(:,r) = ie_dat{s}.rois{r}.dat;
   end
   ie_rois(a:a+size(sub_rois,1)-1,1:size(sub_rois,2))=sub_rois;
   a=a+size(sub_rois,1);
end
dat = ie_rois;
clear  subject_grouping;
% subj group [subj x time]
time=size(ie_dat{1}.rois{1}.dat,1);
% subject_grouping = [ones(r,size(ie_dat,2)*time) (ones(r,size(se_dat,2)*time)*2) (ones(r,size(ne_dat,2)*time)*3)]';
a=1;
for s = 1:size(ie_dat,2)
    subject_grouping(a:a+time-1,1) = [ones(1,time)*s]';
    a=a+time;
end
ieOUT = canlab_connectivity_predict(dat, subject_grouping, 'partialr'); % could use outcome for pattern expression?
% save IE_MultiConnect_Out ieOUT
t =  ieOUT.stats.fdr_thresholded_tvalues;
[stats, handles] = canlab_force_directed_graph(t,'names', names,'rset', ieOUT.parcelindx,'cl',cl)
saveas(gcf, 'IE_connect_graph_fixcolor.png');

%%
clear subject_grouping;clear sub_rois;
a=1;
for s = 1:size(se_dat,2)
   for r = 1:size(se_dat{s}.rois,1)
        sub_rois(:,r) = se_dat{s}.rois{r}.dat;
   end
   se_rois(a:a+size(sub_rois,1)-1,1:size(sub_rois,2))=sub_rois;
   a=a+size(sub_rois,1);
end
dat = se_rois;
% subj group [subj x time]
time=210;
a=1;
for s = 1:size(se_dat,2)
    subject_grouping(a:a+time-1,1) = [ones(1,time)*s]';
    a=a+time;
end
seOUT = canlab_connectivity_predict(dat, subject_grouping, 'partialr'); % could use outcome for pattern expression?
% save SE_MultiConnect_Out seOUT
t =  seOUT.stats.fdr_thresholded_tvalues;
[stats, handles] = canlab_force_directed_graph(t,'names', names,'rset', seOUT.parcelindx,'cl',cl)
saveas(gcf, 'SE_connect_graph_fixcolors.png');

%%
a=1;
clear subject_grouping;clear sub_rois;
for s = 1:size(ne_dat,2)
   for r = 1:size(ne_dat{s}.rois,1)
        sub_rois(:,r) = ne_dat{s}.rois{r}.dat;
   end
   ne_rois(a:a+size(sub_rois,1)-1,1:size(sub_rois,2))=sub_rois;
   a=a+size(sub_rois,1);
end
dat = ne_rois;
% subj group [subj x time]
time=210;
a=1;
for s = 1:size(ne_dat,2)
    subject_grouping(a:a+time-1,1) = [ones(1,time)*s]';
    a=a+time;
end
neOUT = canlab_connectivity_predict(dat, subject_grouping, 'partialr'); % could use outcome for pattern expression?
save NE_MultiConnect_Out neOUT
t =  neOUT.stats.fdr_thresholded_tvalues;
[stats, handles] = canlab_force_directed_graph(t,'names', names,'rset', neOUT.parcelindx,'cl',cl)
saveas(gcf, 'NE_connect_graph.png');

%% use canlab_connectivity_predict to predict outcome (group) OR threat expression! (within group)

load IE_MultiConnect_Out
load SE_MultiConnect_Out
load NE_MultiConnect_Out
%% set up to do binary prediction between the groups
% IE v NE
clear rois;clear dat;clear sub_rois;
close all;
% concat the two dats for comparison
dat=[ie_dat ne_dat];
time=size(ie_dat{1}.rois{1}.dat,1);
a=1;
for s = 1:size(dat,2)
   for r = 1:size(dat{s}.rois,1)
        sub_rois(:,r) = dat{s}.rois{r}.dat;
   end
   rois(a:a+size(sub_rois,1)-1,1:size(sub_rois,2))=sub_rois;
   a=a+size(sub_rois,1);
end
clear  subject_grouping;
clear y;
% subj group [subj x time]
a=1;
for s = 1:size(dat,2)
    subject_grouping(a:a+time-1,1) = [ones(1,time)*s]';
    a=a+time;
end
% set up group identity outcome
% IE = 1
% NE - -1
y=[ones(size(ie_dat,2),1);-ones(size(ne_dat,2),1)];
OUT = canlab_connectivity_predict(rois, subject_grouping, 'partialr','outcome', y,'algo', 'cv_svm', 'folds', 1);
save Compare_IE_v_NE OUT
t =  OUT.stats.fdr_thresholded_tvalues;
[stats, handles] = canlab_force_directed_graph(t,'names', names,'rset', OUT.parcelindx,'cl',cl);

ROC = roc_plot(OUT.PREDICT.pairwise_association.dist_from_hyperplane_xval, logical(OUT.PREDICT.pairwise_association.Y > 0));
    
   
%%
% IE v SE
close all;clear y;
clear rois; clear sub_rois;clear  subject_grouping;
% concat the two dats for comparison
dat=[ie_dat se_dat];
time=size(ie_dat{1}.rois{1}.dat,1);
a=1;
for s = 1:size(dat,2)
   for r = 1:size(dat{s}.rois,1)
        sub_rois(:,r) = dat{s}.rois{r}.dat;
   end
   rois(a:a+size(sub_rois,1)-1,1:size(sub_rois,2))=sub_rois;
   a=a+size(sub_rois,1);
end
% subj group [subj x time]
a=1;
for s = 1:size(dat,2)
    subject_grouping(a:a+time-1,1) = [ones(1,time)*s]';
    a=a+time;
end
% set up group identity outcome
% IE = 1
% SE - -1
y=[ones(size(ie_dat,2),1);-ones(size(se_dat,2),1)];
OUT = canlab_connectivity_predict(rois, subject_grouping, 'partialr','outcome', y,'algo', 'cv_svm', 'folds', 1);
save Compare_IE_v_SE OUT
t =  OUT.stats.fdr_thresholded_tvalues;
[stats, handles] = canlab_force_directed_graph(t,'names', names,'rset', OUT.parcelindx,'cl',cl);

%%
% SE v NE
close all;
clear rois; clear sub_rois;clear  subject_grouping;
% concat the two dats for comparison
dat=[se_dat ne_dat];
time=size(se_dat{1}.rois{1}.dat,1);
a=1;
for s = 1:size(dat,2)
   for r = 1:size(dat{s}.rois,1)
        sub_rois(:,r) = dat{s}.rois{r}.dat;
   end
   rois(a:a+size(sub_rois,1)-1,1:size(sub_rois,2))=sub_rois;
   a=a+size(sub_rois,1);
end
clear  subject_grouping;
% subj group [subj x time]
a=1;
for s = 1:size(dat,2)
    subject_grouping(a:a+time-1,1) = [ones(1,time)*s]';
    a=a+time;
end
% set up group identity outcome
% SE = 1
% NE - -1
y=[ones(size(se_dat,2),1);-ones(size(ne_dat,2),1)];
OUT = canlab_connectivity_predict(rois, subject_grouping, 'partialr','outcome', y,'algo', 'cv_svm', 'folds', 1);
save Compare_SE_v_NE OUT
t =  OUT.stats.fdr_thresholded_tvalues;
[stats, handles] = canlab_force_directed_graph(t,'names', names,'rset', OUT.parcelindx,'cl',cl);


%% load the out and make roc plots

load Compare_IE_v_NE
ROC = roc_plot(OUT.PREDICT.pairwise_association.other_output{2}, logical(OUT.PREDICT.pairwise_association.Y > 0));

load Compare_IE_v_SE
ROC = roc_plot(OUT.PREDICT.pairwise_association.other_output{2}, logical(OUT.PREDICT.pairwise_association.Y > 0));

load Compare_SE_v_NE
ROC = roc_plot(OUT.PREDICT.pairwise_association.other_output{2}, logical(OUT.PREDICT.pairwise_association.Y > 0));
