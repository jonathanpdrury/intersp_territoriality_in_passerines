require(geiger)
require(MCMCglmm)
source('Code0_plmm_functions.R')

phylo<-read.nexus("Tree1.nex") #read in phylogeny
it_dataset<-read.csv("Data1.csv") #read in data


#prepare phylogeny for PLMM:
#trim the phylogeny to only include species that are in the it_dataset dataset
t<-drop.tip(phylo,phylo$tip.label[which(!phylo$tip.label%in%unique(c(as.character(it_dataset$sp1.treename),as.character(it_dataset$sp2.treename))))]) #write this to drop tips that aren't in the dataset
t2<-add.outgroup(t) #adds an arbitrary lineage as an outgroup; this forces root to have a node ID, which it wouldn't otherwise--this is important because some species pairs might share a MRCA at the root
t3<-makeNodeLabel(t2) #convert that phylogeny into a pedigree object for PLMM
INphylo<-inverseA(t3,nodes="ALL")
save(INphylo,file="INphylo.RData") #save file for future use

###The following chunk of code arranges random effects for the phylogeny (not run--these variables already appear in Data1)
##define 'animal' (i.e., a string designating the node name connecting two species) and add it directly to dataset. here,code is taken from Drury et al 2018, Syst Biol:
#animal<-vector()
#for(k in 1:dim(it_dataset)[1]){
#	sp1.treename.node<-INphylo$pedigree[match(it_dataset$sp1.treename[k],INphylo$pedigree[,1]),2]
#	sp1.treename.node.list<-c(sp1.treename.node)
#	while(!is.na(sp1.treename.node)){
#		sp1.treename.node<-INphylo$pedigree[match(sp1.treename.node,INphylo$pedigree[,1]),2]
#		sp1.treename.node.list<-c(sp1.treename.node.list,sp1.treename.node)
#		}
#	sp2.treename.node<-INphylo$pedigree[match(it_dataset$sp2.treename[k],INphylo$pedigree[,1]),2]
#	sp2.treename.node.list<-c(sp2.treename.node)
#	while(!is.na(sp2.treename.node)){
#		sp2.treename.node<-INphylo$pedigree[match(sp2.treename.node,INphylo$pedigree[,1]),2]
#		sp2.treename.node.list<-c(sp2.treename.node.list,sp2.treename.node)
#		}
#	animal<-c(animal,intersect(sp1.treename.node.list,sp2.treename.node.list)[1])	##the first element in common of sp1.treename.node.list and sp2.treename.node.list is their mrca	
#	}
#it_dataset$animal<-animal
#
##next, sort species names between the "sp1" and "sp2" random effects (as in Tobias et al. 2014), so that each species appears more or less the same number of times in "species1" and "species2"
##NB this takes several hours to run; can be saved and reloaded 
#lineages<-as.matrix(data.frame(it_dataset$AOU.1,it_dataset$AOU.2))
#lineages.sorted<-random.effect.sorting(lineages,counter.max=1e6)
#save(lineages.sorted,file="lineages.sorted.RData")
#load("lineages.sorted.RData")
##note lineages.sorted$asymmetry shows the balance for each species between appearing first and second in random effects (==0 if species occurs the same amoutn of time as 'sp1' and 'sp2' in random effects)--this is about as good as we can do, and is mostly evenly spread
##lineages.sorted$sorted.data provides the balanced order for each species
#
##this next block translates the AOU codes back into species names as they appear in the tree and stores them in the sorted order
#species1<-vector()
#for(i in 1:length(lineages.sorted$sorted.data[,1])){
#	if(lineages.sorted$sorted.data[i,1]%in%it_dataset$AOU.1){
#		species1<-c(species1,as.character(it_dataset$sp1.treename[which(it_dataset$AOU.1==lineages.sorted$sorted.data[i,1])][1]))
#		} else {
#		species1<-c(species1,as.character(it_dataset$sp2.treename[which(it_dataset$AOU.2==lineages.sorted$sorted.data[i,1])][1]))
#		}
#	}	
#species2<-vector()
#for(i in 1:length(lineages.sorted$sorted.data[,2])){
#	if(lineages.sorted$sorted.data[i,2]%in%it_dataset$AOU.1){
#		species2<-c(species2,as.character(it_dataset$sp1.treename[which(it_dataset$AOU.1==lineages.sorted$sorted.data[i,2])][1]))
#		} else {
#		species2<-c(species2,as.character(it_dataset$sp2.treename[which(it_dataset$AOU.2==lineages.sorted$sorted.data[i,2])][1]))
#		}
#	}	
#it_dataset$species1<-species1 #match up sorted AOU number to species name
#it_dataset$species2<-species2 #match up sorted AOU number to species name
#
#
#it_dataset$mass.diff.sqrt<-sqrt(abs(log(it_dataset$sp1.mass)-log(it_dataset$sp2.mass)))
#it_dataset$bill.diff.sqrt<-sqrt(abs(log(it_dataset$sp1.male.bill.len)-log(it_dataset$sp2.male.bill.len)))
#it_dataset$same.family<-ifelse(as.character(it_dataset$sp1.family)== as.character(it_dataset$sp2.family),1,0)
#it_dataset$plumage.diff.sq<-it_dataset$Plumage.mean.score^2
#it_dataset$sympatry.sq<-it_dataset$BL.BBS.eBird_Sympatry^2
#it_dataset$simple.habitat<-ifelse(it_dataset$max.habitat.complexity==1,1,0)
#write.csv(it_dataset,file="Data1.csv")

################################################
### Analyses can start here if the above block has been run
################################################
require(MCMCglmm)

#load in data and phylogenetic random effect object
it_dataset<-read.csv("Data1.csv")
load('INphylo.RData')

#For further details:
#see Drury et al. 2018 Syst. Biol., general approach from http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-2-multiple-measurements-model-mcmcglmm
#see also: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q3/004481.html
#and: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q2/016385.html


### INTERSPECIFIC TERRITORIALITY ANALYSES (ALL COMPARISONS)

#Fitting models assessed in Extended Data Table 2

#First, fitting a model for all species pairs with main effects only
prior.m3<-list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)),B=list(mu=c(rep(0,14)),V=diag(14)*(1+pi^2/3)),R=list(V=1,fix=1))
#see justification for this prior in the MCMCglmm course notes section 2.6 (p. 55-56)

m3_chain1<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,hybrid.ok!="NA" &logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m3,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
m3_chain2<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,hybrid.ok!="NA" &logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m3,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
m3_chain3<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,hybrid.ok!="NA" &logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m3,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
m3_chain4<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,hybrid.ok!="NA" &logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m3,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)

##check chain convergence (statistic should be near '1')
m3_combined<-mcmc.list(m3_chain1$Sol,m3_chain2$Sol,m3_chain3$Sol,m3_chain4$Sol)
gelman.diag(m3_combined) #assess chain convergence
plot(m3_combined)  #visualise mixing and convergence in MCMC chains

##calculate phylogenetic signal (i.e., heritability, phylogenetic intra-class correlation, lambda--all are synonyms)
##see: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q3/019060.html

lambda_m31<-m3_chain1$VCV[,'animal']/(rowSums(m3_chain1$VCV)+(pi^2/3))
lambda_m32<-m3_chain2$VCV[,'animal']/(rowSums(m3_chain2$VCV)+(pi^2/3))
lambda_m33<-m3_chain3$VCV[,'animal']/(rowSums(m3_chain3$VCV)+(pi^2/3))
lambda_m34<-m3_chain4$VCV[,'animal']/(rowSums(m3_chain4$VCV)+(pi^2/3))
lambda_m3<-c(as.vector(lambda_m31),as.vector(lambda_m32),as.vector(lambda_m33),as.vector(lambda_m34))

mean(lambda_m3)   #mean lambda (phylogenetic signal) value
quantile(lambda_m3,c(0.025,0.975)) #lambda 95% credibility interval
mean(c(m3_chain1$DIC,m3_chain2$DIC,m3_chain3$DIC,m3_chain4$DIC)) #mean DIC value


#Now, fitting models with interactions between several fixed effects and syntopy:

prior.mx0<-list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)),B=list(mu=c(rep(0,17)),V=diag(17)*(1+pi^2/3)),R=list(V=1,fix=1))

#model with all species pairs and both syntopy*hybrid and syntopy*mass interactions
mx0_chain1<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt)+ scale(logBBSe.syntopyOu)*scale(proportion.shared.axes) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx0,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
mx0_chain2<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt)+ scale(logBBSe.syntopyOu)*scale(proportion.shared.axes) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx0,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
mx0_chain3<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt)+ scale(logBBSe.syntopyOu)*scale(proportion.shared.axes) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx0,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
mx0_chain4<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt)+ scale(logBBSe.syntopyOu)*scale(proportion.shared.axes) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx0,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)

mx0_combined<-mcmc.list(mx0_chain1$Sol,mx0_chain2$Sol,mx0_chain3$Sol,mx0_chain4$Sol)
gelman.diag(mx0_combined) #assess chain convergence
plot(mx0_combined) #visualise mixing and convergence in MCMC chains
mean(c(mx0_chain1$DIC,mx0_chain2$DIC,mx0_chain3$DIC,mx0_chain4$DIC)) #mean DIC value

prior.mx0B<-list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)),B=list(mu=c(rep(0,18)),V=diag(18)*(1+pi^2/3)),R=list(V=1,fix=1))

#model with all species pairs and both syntopy*hybrid and syntopy*mass interactions
mx0B_chain1<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt)+ scale(logBBSe.syntopyOu)*scale(proportion.shared.axes) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester)+scale(logBBSe.syntopyOu)*factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx0B,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
mx0B_chain2<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt)+ scale(logBBSe.syntopyOu)*scale(proportion.shared.axes) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester)+scale(logBBSe.syntopyOu)*factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx0B,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
mx0B_chain3<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt)+ scale(logBBSe.syntopyOu)*scale(proportion.shared.axes) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester)+scale(logBBSe.syntopyOu)*factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx0B,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
mx0B_chain4<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt)+ scale(logBBSe.syntopyOu)*scale(proportion.shared.axes) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester)+scale(logBBSe.syntopyOu)*factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx0B,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)

mx0B_combined<-mcmc.list(mx0B_chain1$Sol,mx0B_chain2$Sol,mx0B_chain3$Sol,mx0B_chain4$Sol)
gelman.diag(mx0B_combined) #assess chain convergence
plot(mx0B_combined) #visualise mixing and convergence in MCMC chains
mean(c(mx0B_chain1$DIC,mx0B_chain2$DIC,mx0B_chain3$DIC,mx0B_chain4$DIC)) #mean DIC value


prior.mx1<-list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)),B=list(mu=c(rep(0,16)),V=diag(16)*(1+pi^2/3)),R=list(V=1,fix=1))

#model with all species pairs and both syntopy*hybrid and syntopy*mass interactions
mx1_chain1<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx1,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
mx1_chain2<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx1,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
mx1_chain3<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx1,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
mx1_chain4<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(logBBSe.syntopyOu)*scale(mass.diff.sqrt) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx1,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)

mx1_combined<-mcmc.list(mx1_chain1$Sol,mx1_chain2$Sol,mx1_chain3$Sol,mx1_chain4$Sol)
gelman.diag(mx1_combined) #assess chain convergence
plot(mx1_combined) #visualise mixing and convergence in MCMC chains
mean(c(mx1_chain1$DIC,mx1_chain2$DIC,mx1_chain3$DIC,mx1_chain4$DIC)) #mean DIC value

prior.mx2<-list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)),B=list(mu=c(rep(0,15)),V=diag(15)*(1+pi^2/3)),R=list(V=1,fix=1))

#model with all species pairs and just the syntopy*hybrid interaction
#(Presented in Extended Data Table 3)
mx2_chain1<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx2,nitt=2000000,burnin=20000,thin=1000,verbose=FALSE,slice=TRUE)
mx2_chain1<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx2,nitt=2000000,burnin=20000,thin=1000,verbose=FALSE,slice=TRUE)
mx2_chain1<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx2,nitt=2000000,burnin=20000,thin=1000,verbose=FALSE,slice=TRUE)
mx2_chain1<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu) + scale(logBBSe.syntopyOu)*factor(hybrid.ok) + scale(sympatry.sq) + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset,logBBSe.syntopyOu!="NA"), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.mx2,nitt=2000000,burnin=20000,thin=1000,verbose=FALSE,slice=TRUE)

#check chain convergence (statistic should be near '1')
mx2_combined<-mcmc.list(mx2_chain1$Sol,mx2_chain2$Sol,mx2_chain3$Sol,mx2_chain4$Sol)
gelman.diag(mx2_combined) #assess chain convergence
plot(mx2_combined) #visualise mixing and convergence in MCMC chains

lambda1<-mx2_chain1$VCV[,'animal']/(rowSums(mx2_chain1$VCV)+(pi^2/3))
lambda2<-mx2_chain2$VCV[,'animal']/(rowSums(mx2_chain2$VCV)+(pi^2/3))
lambda3<-mx2_chain3$VCV[,'animal']/(rowSums(mx2_chain3$VCV)+(pi^2/3))
lambda4<-mx2_chain4$VCV[,'animal']/(rowSums(mx2_chain4$VCV)+(pi^2/3))
lambda_mx2<-c(as.vector(lambda1),as.vector(lambda2),as.vector(lambda3),as.vector(lambda4))

mean(lambda_mx2)  #mean lambda (phylogenetic signal) value
quantile(lambda_mx2,c(0.025,0.975)) #lambda 95% credibility interval
mean(c(mx2_chain1$DIC,mx2_chain2$DIC,mx2_chain3$DIC,mx2_chain4$DIC)) #mean DIC value



### INTRAFAMILY INTERSPECIFIC TERRITORIALITY ANALYSES
## Analyses in Extended Data Table 4a:

prior.m4a<-list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)),B=list(mu=c(rep(0,14)),V=diag(14)*(1+pi^2/3)),R=list(V=1,fix=1))
#see justification for this prior in the course notes section 2.6 (p. 55-56)

m4a_chain1<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu)+ scale(sympatry.sq)  + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset, logBBSe.syntopyOu!="NA" & same.family==1), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m4a,nitt=2000000,burnin=20000,thin=1000,verbose=TRUE,slice=TRUE)
m4a_chain2<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu)+ scale(sympatry.sq)  + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset, logBBSe.syntopyOu!="NA" & same.family==1), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m4a,nitt=2000000,burnin=20000,thin=1000,verbose=FALSE,slice=TRUE)
m4a_chain3<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu)+ scale(sympatry.sq)  + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset, logBBSe.syntopyOu!="NA" & same.family==1), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m4a,nitt=2000000,burnin=20000,thin=1000,verbose=FALSE,slice=TRUE)
m4a_chain4<-MCMCglmm(factor(IT.classif) ~ factor(hybrid.ok) + scale(logBBSe.syntopyOu)+ scale(sympatry.sq)  + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset, logBBSe.syntopyOu!="NA" & same.family==1), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m4a,nitt=2000000,burnin=20000,thin=1000,verbose=FALSE,slice=TRUE)

m4a_combined<-mcmc.list(m4a_chain1$Sol,m4a_chain2$Sol,m4a_chain3$Sol,m4a_chain4$Sol)
gelman.diag(m4a_combined) #assess chain convergence
plot(m4a_combined) #visualise mixing and convergence in MCMC chains

lambda_m4a1<-m4a_chain1$VCV[,'animal']/(rowSums(m4a_chain1$VCV)+(pi^2/3))
lambda_m4a2<-m4a_chain2$VCV[,'animal']/(rowSums(m4a_chain2$VCV)+(pi^2/3))
lambda_m4a3<-m4a_chain3$VCV[,'animal']/(rowSums(m4a_chain3$VCV)+(pi^2/3))
lambda_m4a4<-m4a_chain4$VCV[,'animal']/(rowSums(m4a_chain4$VCV)+(pi^2/3))
lambda_m4a<-c(as.vector(lambda1),as.vector(lambda2),as.vector(lambda3),as.vector(lambda4))

mean(lambda_m4a) #mean lambda (phylogenetic signal) value
quantile(lambda_m4a,c(0.025,0.975)) #lambda 95% credibility interval
mean(c(m4a_chain1$DIC,m4a_chain2$DIC,m4a_chain3$DIC,m4a_chain4$DIC)) #mean DIC value

### INTERFAMILY INTERSPECIFIC TERRITORIALITY ANALYSES
## Analyses in Extended Data Table 4b:

prior.m4b<-list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)),B=list(mu=c(rep(0,13)),V=diag(13)*(1+pi^2/3)),R=list(V=1,fix=1))
#see justification for this prior in the course notes section 2.6 (p. 55-56)

m4b_chain1<-MCMCglmm(factor(IT.classif) ~ scale(logBBSe.syntopyOu)+ scale(sympatry.sq)  + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset, logBBSe.syntopyOu!="NA" & same.family==0), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m4b,nitt=10000000,burnin=50000,thin=1000,verbose=TRUE,slice=TRUE)
m4b_chain2<-MCMCglmm(factor(IT.classif) ~ scale(logBBSe.syntopyOu)+ scale(sympatry.sq)  + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset, logBBSe.syntopyOu!="NA" & same.family==0), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m4b,nitt=10000000,burnin=50000,thin=1000,verbose=TRUE,slice=TRUE)
m4b_chain3<-MCMCglmm(factor(IT.classif) ~ scale(logBBSe.syntopyOu)+ scale(sympatry.sq)  + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset, logBBSe.syntopyOu!="NA" & same.family==0), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m4b,nitt=10000000,burnin=50000,thin=1000,verbose=TRUE,slice=TRUE)
m4b_chain4<-MCMCglmm(factor(IT.classif) ~ scale(logBBSe.syntopyOu)+ scale(sympatry.sq)  + scale(plumage.diff.sq) + scale(patristic.distance) + scale(proportion.shared.axes)+scale(mass.diff.sqrt)+scale(bill.diff.sqrt)+factor(simple.habitat)+scale(song.pcdists)+scale(mean.max.xcor)+factor(same.territory.type)+factor(both.hole.nester), random = ~animal + species1 + species2, data=subset(it_dataset, logBBSe.syntopyOu!="NA" & same.family==0), family = "categorical", ginverse = list(animal = INphylo$Ainv),prior=prior.m4b,nitt=10000000,burnin=50000,thin=1000,verbose=TRUE,slice=TRUE)

m4b_combined<-mcmc.list(m4b_chain1$Sol,m4b_chain2$Sol,m4b_chain3$Sol,m4b_chain4$Sol)
gelman.diag(m4b_combined) #assess chain convergence
plot(m4b_combined) #visualise mixing and convergence in MCMC chains

lambda_m4b1<-m4b_chain1$VCV[,'animal']/(rowSums(m4b_chain1$VCV)+(pi^2/3))
lambda_m4b2<-m4b_chain2$VCV[,'animal']/(rowSums(m4b_chain2$VCV)+(pi^2/3))
lambda_m4b3<-m4b_chain3$VCV[,'animal']/(rowSums(m4b_chain3$VCV)+(pi^2/3))
lambda_m4b4<-m4b_chain4$VCV[,'animal']/(rowSums(m4b_chain4$VCV)+(pi^2/3))
lambda_m4b<-c(as.vector(lambda1),as.vector(lambda2),as.vector(lambda3),as.vector(lambda4))

mean(lambda_m4b) #mean lambda (phylogenetic signal) value
quantile(lambda_m4b,c(0.025,0.975)) #lambda 95% credibility interval
mean(c(m4b_chain1$DIC,m4b_chain2$DIC,m4b_chain3$DIC,m4b_chain4$DIC)) #mean DIC value

