%% Attenuating neural threat expression with imagination
%  Reddan, Wager, & Schiller, Neuron 2018
% ------------------------------------------------------------
% this script uses code from
% https://github.com/canlab
% SPM8
% ------------------------------------------------------------
% Figure 5. Imagined and real extinction recruit the vmPFC. A. vmPFC activation 
% across experimental phases. Left. Depiction of vmPFC anatomical mask used in 
% this analysis. Due to the role of the vmPFC during extinction, the vmPFC was 
% an a priori region of interest. This bilateral mask was created with Neurosynth. 
% Right. Average differential (CS+ > CS-) BOLD signal in the vmPFC rose over the 
% course of extinction, peaked in the third bin, and then decreased in the imagined 
% and real extinction groups. The no extinction group exhibited little to no 
% change in activation across binned trials.  The spatial voxel patterns during 
% the third extinction time point are displayed, revealing (nonsignificant) 
% similarity in the activity patterns between the imagined and real groups. 
% During late re-extinction, vmPFC activity increased in the imagined and real 
% groups, but not in the no extinction group, which demonstrated threat responses 
% during this phase. Error bars reflect standard error of the mean.

%% A - Univariate vmPFC activation across each phase

load IE_Figure5A_dat
t = table(group,ALL(:,1),ALL(:,2),ALL(:,3),ALL(:,4),ALL(:,5),ALL(:,6),ALL(:,7),...
'VariableNames',{'group','acq','ext1','ext2','ext3','ext4','reext1','reext2'});
Meas = dataset([1 2 3 4 5 6 7]','VarNames',{'Measurements'});
rm = fitrm(t,'acq-reext2~group','WithinDesign',Meas);
ranovatbl = ranova(rm)

NE=ALL(1:24,1:end);
SE=ALL(25:46,1:end);
IE=ALL(47:end,1:end);
Plot_Title='vmPFC activation across all Phases';
Y_Label=sprintf('Differential (CS+ > CS-) \n Beta Weights');
ie_get_lineplot(IE, SE, NE, Plot_Title, Y_Label)