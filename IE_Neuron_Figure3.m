%% Attenuating neural threat expression with imagination
%  Reddan, Wager, & Schiller, Neuron 2018
% ------------------------------------------------------------
% this script uses code from
% https://github.com/canlab
% SPM8
% ------------------------------------------------------------
% Figure 3. Imagined Extinction Reduces Neural and Physiological Threat Expression. 
% A. Neural Threat-Predictive Pattern Expression during the Recovery Test. The threat-
% predictive pattern was applied to brain activity during the late re-extinction 
% recovery test. Only the no extinction group (?none?, N=24) demonstrated threat 
% recovery during late re-extinction, indicated by a comparison of average threat-
% predictive pattern expression values between groups. This may indicate that both 
% imagined (N=20) and real (N=22) extinction successfully generated an extinction memory. 
% Error bars reflect standard error of the mean. B. SCRs during the Recovery Test. 
% Physiological findings, indicated by SCR, complimented the neural findings. Relative 
% to no extinction (N=13), imagined (N=12) and real extinction (N=17) showed a reduction 
% in threat-related physiological arousal during the recovery test. Error bars reflect 
% standard error of the mean. Participants (N=24) who did not demonstrate a discriminatory 
% SCR during acquisition, defined as greater SCR to the CS+ relative to the CS- during 
% either the first or last half of threat-acquisition on average, were excluded from SCR 
% analysis because we were unable to determine if their SCR were representative of threat-
% related arousal. C. Correlation between SCRs and Threat-Predictive Pattern Expression 
% during the Recovery Test. SCRs were positively correlated with expression of the threat-
% predictive pattern (rho = 0.34, p = 0.03), indicating a link between the neural and 
% physiological expressions of threat. 

%% A - Neural threat-predictive pattern expression during recovery test
load IE_Figure3A_dat;
Plot_Title='Whole Brain Threat Pattern Expression during Recovery Test';
Y_Label='Average Pattern Expression';
ie_getstats(IE, SE, NE, Plot_Title, Y_Label);
% -------------------- output --------------------
% ---- Whole Brain Threat Pattern Expression during Recovery Test----
% sig anova
% F(2,63)=3.535460, p=0.035058, Eta Sq=0.100911, means IE: -0.202705 SE: -0.114120 NE: 0.548541
% IE v SE: t(31.842222)=-0.273660,p=0.786113, CI=[-0.748075 0.570906], STD=[1.234971 0.792197], Hedges g=-0.084661
% SE v NE: t(42.248941)=-2.404495,p=0.020655, CI=[-1.218731 -0.106590], STD=[0.792197 1.066864], Hedges g=-0.688641
% IE v NE: t(37.878622)=-2.136129,p=0.039193, CI=[-1.463270 -0.039220], STD=[1.234971 1.066864], Hedges g=-0.643777

% Supplementary Figure 2 A
ie_getstats_violin(IE, SE, NE, Plot_Title, Y_Label);

%% B - SCR during recovery test
load IE_Figure3B_dat;
Plot_Title='SCR during Recovery Test';
Y_Label='microsemens';
ie_getstats(IE, SE, NE, Plot_Title, Y_Label);
% -------------------- output --------------------
% ---- SCR during Recovery Test----
% sig anova
% F(2,39)=3.609121, p=0.036467, Eta Sq=0.156177, means IE: -0.000381 SE: -0.004815 NE: 0.034581
% IE v SE: t(26.468711)=0.317683,p=0.753218, CI=[-0.024229 0.033096], STD=[0.028098 0.046825], Hedges g=0.107032
% SE v NE: t(26.206415)=-2.306155,p=0.029263, CI=[-0.074497 -0.004295], STD=[0.046825 0.046012], Hedges g=-0.824712
% IE v NE: t(20.079422)=-2.312186,p=0.031495, CI=[-0.066496 -0.003429], STD=[0.028098 0.046012], Hedges g=-0.878209

% Supplementary Figure 2 B
ie_getstats_violin(IE, SE, NE, Plot_Title, Y_Label);

%% C - Correlation between SCR & Neural threat-predictive pattern expression during recovery test

load IE_Figure3C_dat;
[r p]=corr(ALL42, 'Type','Spearman'); % r = 0.3362, p=0.03
x=ALL42(:,1);%scr
y=ALL42(:,2);%PE
bls = regress(y,[ones(length(x),1) x]);
scatter(x,y,'filled'); grid on; hold on
plot(x,bls(1)+bls(2)*x,'r','LineWidth',2);

