require(MCMCglmm)

#load in data and phylogenetic random effect object

it_dataset<-read.csv("Data1.csv")
load('INphylo.RData') #see Code1_interspecific_territoriality_PLMMs.R for instructions on constructing this file

it_dataset$hybrid.cavity<-apply(cbind(it_dataset$hybrid.ok,it_dataset$both.hole.nester),1,max)

#define the prior (a generic, non-informative prior):
prior<-list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)),R=list(V=1,nu=0.02))

### INTRAFAMILY PLUMAGE ANALYSES

## Analyses in Table S9:

#first, a model with main effects and no interactions for all intrafamily IT species pairs
m.p_chain1<-MCMCglmm(scale(plumage.diff.sq) ~ factor(hybrid.ok) + factor(IT.classif)  + scale(logBBSe.syntopyOu)+ scale(proportion.shared.axes) + scale(mass.diff.sqrt)+ factor(both.hole.nester)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA"), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p_chain2<-MCMCglmm(scale(plumage.diff.sq) ~ factor(hybrid.ok) + factor(IT.classif)  + scale(logBBSe.syntopyOu)+ scale(proportion.shared.axes) + scale(mass.diff.sqrt)+ factor(both.hole.nester)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA"), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p_chain3<-MCMCglmm(scale(plumage.diff.sq) ~ factor(hybrid.ok) + factor(IT.classif) + scale(logBBSe.syntopyOu) + scale(proportion.shared.axes) + scale(mass.diff.sqrt)+ factor(both.hole.nester)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA"), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p_chain4<-MCMCglmm(scale(plumage.diff.sq) ~ factor(hybrid.ok) + factor(IT.classif)  + scale(logBBSe.syntopyOu)+ scale(proportion.shared.axes) + scale(mass.diff.sqrt)+ factor(both.hole.nester)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA"), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)

m.p_combined<-mcmc.list(m.p_chain1$Sol,m.p_chain2$Sol,m.p_chain3$Sol,m.p_chain4$Sol)
gelman.diag(m.p_combined) #assess chain convergence
plot(m.p_combined) #visualise mixing and convergence in MCMC chains

lambda_m.p1<-m.p_chain1$VCV[,'animal']/(rowSums(m.p_chain1$VCV)+(pi^2/3))
lambda_m.p2<-m.p_chain2$VCV[,'animal']/(rowSums(m.p_chain2$VCV)+(pi^2/3))
lambda_m.p3<-m.p_chain3$VCV[,'animal']/(rowSums(m.p_chain3$VCV)+(pi^2/3))
lambda_m.p4<-m.p_chain4$VCV[,'animal']/(rowSums(m.p_chain4$VCV)+(pi^2/3))
lambda_m.p<-c(as.vector(lambda_m.p1),as.vector(lambda_m.p2),as.vector(lambda_m.p3),as.vector(lambda_m.p4))

mean(lambda_m.p) #mean lambda (phylogenetic signal) value
quantile(lambda_m.p,c(0.025,0.975)) #lambda 95% credibility interval
mean(c(m.p_chain1$DIC,m.p_chain2$DIC,m.p_chain3$DIC,m.p_chain4$DIC)) #mean DIC value

#second, a model with main effects and  interactions for all intrafamily IT pairs
m.p.x_chain1<-MCMCglmm(scale(plumage.diff.sq) ~ factor(hybrid.ok)*factor(IT.classif) + scale(logBBSe.syntopyOu)*factor(IT.classif) + scale(proportion.shared.axes)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif)+ factor(both.hole.nester)*factor(IT.classif)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA"), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p.x_chain2<-MCMCglmm(scale(plumage.diff.sq) ~ factor(hybrid.ok)*factor(IT.classif) + scale(logBBSe.syntopyOu)*factor(IT.classif) + scale(proportion.shared.axes)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif)+ factor(both.hole.nester)*factor(IT.classif)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA"), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p.x_chain3<-MCMCglmm(scale(plumage.diff.sq) ~ factor(hybrid.ok)*factor(IT.classif) + scale(logBBSe.syntopyOu)*factor(IT.classif) + scale(proportion.shared.axes)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif)+ factor(both.hole.nester)*factor(IT.classif)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA"), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p.x_chain4<-MCMCglmm(scale(plumage.diff.sq) ~ factor(hybrid.ok)*factor(IT.classif) + scale(logBBSe.syntopyOu)*factor(IT.classif) + scale(proportion.shared.axes)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif)+ factor(both.hole.nester)*factor(IT.classif)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA"), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)

m.p.x_combined<-mcmc.list(m.p.x_chain1$Sol,m.p.x_chain2$Sol,m.p.x_chain3$Sol,m.p.x_chain4$Sol)
gelman.diag(m.p.x_combined) #assess chain convergence
plot(m.p.x_combined) #visualise mixing and convergence in MCMC chains

lambda_m.p.x1<-m.p.x_chain1$VCV[,'animal']/(rowSums(m.p.x_chain1$VCV)+(pi^2/3))
lambda_m.p.x2<-m.p.x_chain2$VCV[,'animal']/(rowSums(m.p.x_chain2$VCV)+(pi^2/3))
lambda_m.p.x3<-m.p.x_chain3$VCV[,'animal']/(rowSums(m.p.x_chain3$VCV)+(pi^2/3))
lambda_m.p.x4<-m.p.x_chain4$VCV[,'animal']/(rowSums(m.p.x_chain4$VCV)+(pi^2/3))
lambda_m.p.x<-c(as.vector(lambda_m.p.x1),as.vector(lambda_m.p.x2),as.vector(lambda_m.p.x3),as.vector(lambda_m.p.x4))

mean(lambda_m.p.x) #mean lambda (phylogenetic signal) value
quantile(lambda_m.p.x,c(0.025,0.975)) #lambda 95% credibility interval
mean(c(m.p.x_chain1$DIC,m.p.x_chain2$DIC,m.p.x_chain3$DIC,m.p.x_chain4$DIC)) #mean DIC value


## Analyses in Table S10:

#A model with main effects and no interactions, just non-hybridizing, non-cavity nesting intrafamily IT species pairs
m.p.a_chain1<-MCMCglmm(scale(plumage.diff.sq) ~  scale(patristic.distance)+scale(logBBSe.syntopyOu)+factor(IT.classif)+ scale(proportion.shared.axes) + scale(mass.diff.sqrt), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA" & hybrid.cavity==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p.a_chain2<-MCMCglmm(scale(plumage.diff.sq) ~  scale(patristic.distance)+scale(logBBSe.syntopyOu)+factor(IT.classif)+ scale(proportion.shared.axes) + scale(mass.diff.sqrt), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA" & hybrid.cavity==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p.a_chain3<-MCMCglmm(scale(plumage.diff.sq) ~  scale(patristic.distance)+scale(logBBSe.syntopyOu)+factor(IT.classif)+ scale(proportion.shared.axes) + scale(mass.diff.sqrt), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA" & hybrid.cavity==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p.a_chain4<-MCMCglmm(scale(plumage.diff.sq) ~  scale(patristic.distance)+scale(logBBSe.syntopyOu)+factor(IT.classif)+ scale(proportion.shared.axes) + scale(mass.diff.sqrt), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA" & hybrid.cavity==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)

m.p.a_combined<-mcmc.list(m.p.a_chain1$Sol,m.p.a_chain2$Sol,m.p.a_chain3$Sol,m.p.a_chain4$Sol)
gelman.diag(m.p.a_combined)#assess chain convergence
plot(m.p.a_combined) #visualise mixing and convergence in MCMC chains

lambda_m.p.a1<-m.p.a_chain1$VCV[,'animal']/(rowSums(m.p.a_chain1$VCV)+(pi^2/3))
lambda_m.p.a2<-m.p.a_chain2$VCV[,'animal']/(rowSums(m.p.a_chain2$VCV)+(pi^2/3))
lambda_m.p.a3<-m.p.a_chain3$VCV[,'animal']/(rowSums(m.p.a_chain3$VCV)+(pi^2/3))
lambda_m.p.a4<-m.p.a_chain4$VCV[,'animal']/(rowSums(m.p.a_chain4$VCV)+(pi^2/3))
lambda_m.p.a<-c(as.vector(lambda_m.p.a1),as.vector(lambda_m.p.a2),as.vector(lambda_m.p.a3),as.vector(lambda_m.p.a4))

mean(lambda_m.p.a) #mean lambda (phylogenetic signal) value
quantile(lambda_m.p.a,c(0.025,0.975)) #lambda 95% credibility interval
range(c(m.p.a_chain1$DIC,m.p.a_chain2$DIC,m.p.a_chain3$DIC,m.p.a_chain4$DIC)) #mean DIC value


#Now, a model with main effects and also interactions, for just non-hybridizing, non-cavity nesting intrafamily IT species pairs

m.p.a.x_chain1<-MCMCglmm(scale(plumage.diff.sq) ~  scale(patristic.distance)+scale(logBBSe.syntopyOu)*factor(IT.classif)+ scale(proportion.shared.axes)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA" & hybrid.cavity==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p.a.x_chain2<-MCMCglmm(scale(plumage.diff.sq) ~  scale(patristic.distance)+scale(logBBSe.syntopyOu)*factor(IT.classif)+ scale(proportion.shared.axes)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA" & hybrid.cavity==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p.a.x_chain3<-MCMCglmm(scale(plumage.diff.sq) ~  scale(patristic.distance)+scale(logBBSe.syntopyOu)*factor(IT.classif)+ scale(proportion.shared.axes)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA" & hybrid.cavity==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.p.a.x_chain4<-MCMCglmm(scale(plumage.diff.sq) ~  scale(patristic.distance)+scale(logBBSe.syntopyOu)*factor(IT.classif)+ scale(proportion.shared.axes)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==1 & logBBSe.syntopyOu!="NA" & hybrid.cavity==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)

m.p.a.x_combined<-mcmc.list(m.p.a.x_chain1$Sol,m.p.a.x_chain2$Sol,m.p.a.x_chain3$Sol,m.p.a.x_chain4$Sol)
gelman.diag(m.p.a.x_combined) #assess chain convergence
plot(m.p.a.x_combined) #visualise mixing and convergence in MCMC chains

lambda_m.p.a.x1<-m.p.a.x_chain1$VCV[,'animal']/(rowSums(m.p.a.x_chain1$VCV)+(pi^2/3))
lambda_m.p.a.x2<-m.p.a.x_chain2$VCV[,'animal']/(rowSums(m.p.a.x_chain2$VCV)+(pi^2/3))
lambda_m.p.a.x3<-m.p.a.x_chain3$VCV[,'animal']/(rowSums(m.p.a.x_chain3$VCV)+(pi^2/3))
lambda_m.p.a.x4<-m.p.a.x_chain4$VCV[,'animal']/(rowSums(m.p.a.x_chain4$VCV)+(pi^2/3))
lambda_m.p.a.x<-c(as.vector(lambda_m.p.a.x1),as.vector(lambda_m.p.a.x2),as.vector(lambda_m.p.a.x3),as.vector(lambda_m.p.a.x4))

mean(lambda_m.p.a.x) #mean lambda (phylogenetic signal) value
quantile(lambda_m.p.a.x,c(0.025,0.975)) #lambda 95% credibility interval
mean(c(m.p.a.x_chain1$DIC,m.p.a.x_chain2$DIC,m.p.a.x_chain3$DIC,m.p.a.x_chain4$DIC))#mean DIC value


### INTERFAMILY SONG ANALYSES

## Analyses in Table S11:

##Fit a model with main effects and no interactions for all interfamily IT species pairs:

m.s_chain1<-MCMCglmm(scale(song.pcdists.std) ~  scale(proportion.shared.axes) + factor(IT.classif) + scale(logBBSe.syntopyOu)+ scale(mass.diff.sqrt)+ factor(both.hole.nester)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.s_chain2<-MCMCglmm(scale(song.pcdists.std) ~  scale(proportion.shared.axes) + factor(IT.classif)+ scale(logBBSe.syntopyOu) + scale(mass.diff.sqrt)+ factor(both.hole.nester)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.s_chain3<-MCMCglmm(scale(song.pcdists.std) ~  scale(proportion.shared.axes) + factor(IT.classif) + scale(logBBSe.syntopyOu)+ scale(mass.diff.sqrt)+ factor(both.hole.nester)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.s_chain4<-MCMCglmm(scale(song.pcdists.std) ~  scale(proportion.shared.axes) + factor(IT.classif) + scale(logBBSe.syntopyOu)+ scale(mass.diff.sqrt)+ factor(both.hole.nester)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)

m.s_combined<-mcmc.list(m.s_chain1$Sol,m.s_chain2$Sol,m.s_chain3$Sol,m.s_chain4$Sol)
gelman.diag(m.s_combined) #assess chain convergence
plot(m.s_combined) #visualise mixing and convergence in MCMC chains

lambda_m.s1<-m.s_chain1$VCV[,'animal']/(rowSums(m.s_chain1$VCV)+(pi^2/3))
lambda_m.s2<-m.s_chain2$VCV[,'animal']/(rowSums(m.s_chain2$VCV)+(pi^2/3))
lambda_m.s3<-m.s_chain3$VCV[,'animal']/(rowSums(m.s_chain3$VCV)+(pi^2/3))
lambda_m.s4<-m.s_chain4$VCV[,'animal']/(rowSums(m.s_chain4$VCV)+(pi^2/3))
lambda_m.s<-c(as.vector(lambda_m.s1),as.vector(lambda_m.s2),as.vector(lambda_m.s3),as.vector(lambda_m.s4))

mean(lambda_m.s) #mean lambda (phylogenetic signal) value
quantile(lambda_m.s,c(0.025,0.975)) #lambda 95% credibility interval
mean(c(m.s_chain1$DIC,m.s_chain2$DIC,m.s_chain3$DIC,m.s_chain4$DIC)) #mean DIC value


m.s.x_chain1<-MCMCglmm(scale(song.pcdists.std) ~  scale(proportion.shared.axes)*factor(IT.classif)+ scale(logBBSe.syntopyOu)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif)+ factor(both.hole.nester)*factor(IT.classif)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.s.x_chain2<-MCMCglmm(scale(song.pcdists.std) ~  scale(proportion.shared.axes)*factor(IT.classif)+ scale(logBBSe.syntopyOu)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif)+ factor(both.hole.nester)*factor(IT.classif)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.s.x_chain3<-MCMCglmm(scale(song.pcdists.std) ~  scale(proportion.shared.axes)*factor(IT.classif)+ scale(logBBSe.syntopyOu)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif)+ factor(both.hole.nester)*factor(IT.classif)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)
m.s.x_chain4<-MCMCglmm(scale(song.pcdists.std) ~  scale(proportion.shared.axes)*factor(IT.classif)+ scale(logBBSe.syntopyOu)*factor(IT.classif) + scale(mass.diff.sqrt)*factor(IT.classif)+ factor(both.hole.nester)*factor(IT.classif)+ scale(patristic.distance), random = ~animal + species1 + species2, data=subset(it_dataset,same.family==0), family = "gaussian", ginverse = list(animal = INphylo$Ainv),prior=prior,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE)

m.s.x_combined<-mcmc.list(m.s.x_chain1$Sol,m.s.x_chain2$Sol,m.s.x_chain3$Sol,m.s.x_chain4$Sol)
gelman.diag(m.s.x_combined) #assess chain convergence
plot(m.s.x_combined) #visualise mixing and convergence in MCMC chains

lambda1<-m.s.x_chain1$VCV[,'animal']/(rowSums(m.s.x_chain1$VCV)+(pi^2/3))
lambda2<-m.s.x_chain2$VCV[,'animal']/(rowSums(m.s.x_chain2$VCV)+(pi^2/3))
lambda3<-m.s.x_chain3$VCV[,'animal']/(rowSums(m.s.x_chain3$VCV)+(pi^2/3))
lambda4<-m.s.x_chain4$VCV[,'animal']/(rowSums(m.s.x_chain4$VCV)+(pi^2/3))
lambda<-c(as.vector(lambda1),as.vector(lambda2),as.vector(lambda3),as.vector(lambda4))

mean(lambda) #mean lambda (phylogenetic signal) value
quantile(lambda,c(0.025,0.975)) #lambda 95% credibility interval
mean(c(m.s.x_chain1$DIC,m.s.x_chain2$DIC,m.s.x_chain3$DIC,m.s.x_chain4$DIC)) #mean DIC value



