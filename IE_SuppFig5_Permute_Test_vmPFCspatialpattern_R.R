## Permutation Testing
# http://thomasleeper.com/Rcourse/Tutorials/permutationtests.html

dat=dat_iesene_vmpfc_3rdtp_forpermute
N=nrow(dat)
grp=rep(1:3,c(20,22,24)) #actual treatment vector
# take avg and get correlations btwn groups
ie<-colMeans(dat[grp==1,])
se<-colMeans(dat[grp==2,])
ne<-colMeans(dat[grp==3,])
iese<-cor(ie,se,method = "pearson")
iene<-cor(ie,ne,method = "pearson")
nese<-cor(ne,se,method = "pearson")

# resample group labels
# take group mean
# get correlation
# permute

# resample without replacement and calculate the corr again [or with replacement - bootstrap, uneven N avoid skew]
# repeat this process a large number of times, build permutation distribution
iese_dist <- 0
iene_dist <- 0
nese_dist <- 0
nperms=10000
for(i in 1:nperms){
  s <- sample(grp, length(grp), TRUE) #permuted treatment vector
  s_iese <- cor(colMeans(dat[s==1,]),colMeans(dat[s==2,]),method = "pearson")
  s_iene <- cor(colMeans(dat[s==1,]),colMeans(dat[s==3,]),method = "pearson")
  s_nese <- cor(colMeans(dat[s==3,]),colMeans(dat[s==2,]),method = "pearson")
  iese_dist[i] = s_iese
  iene_dist[i] = s_iene
  nese_dist[i] = s_nese
}
hist(fisherz(iese_dist), xlim = c(-1, 1), col = "black", breaks = 100) # look at distribution using hist
abline(v = fisherz(iese), col = "blue", lwd = 2) # draw a vertical line for our observed difference
sum(abs(fisherz(iese_dist)) > abs(fisherz(iese)))/nperms
# [1] 0.2482

hist(fisherz(iene_dist), xlim = c(-1, 1), col = "black", breaks = 100) # look at distribution using hist
abline(v = fisherz(iene), col = "blue", lwd = 2) # draw a vertical line for our observed difference
sum(abs(fisherz(iene_dist)) > abs(fisherz(iene)))/nperms
# [1] 0.7578

hist(fisherz(nese_dist), xlim = c(-1, 1), col = "black", breaks = 100) # look at distribution using hist
abline(v = fisherz(nese), col = "blue", lwd = 2) # draw a vertical line for our observed difference
sum(abs(fisherz(nese_dist)) > abs(fisherz(nese)))/nperms # two-tailed test
# [1] 0.6609