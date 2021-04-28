# Data analysis pipeline for Soil dataset
library(ape)
library(stringr)
library(magrittr)
library(phylofactor)
### OTU Table
Data <- read.csv('otu_table_wTax_40000_filt2.txt',sep='\t',row.names=1,quote='#',)
taxonomy <- Data[,dim(Data)[2]]
Tax <- data.frame('OTU'=rownames(Data),'taxonomy'=taxonomy)
Data <- Data[,1:(dim(Data)[2]-1)]

otus <- as.list(rownames(Data)) %>% sapply(.,FUN = function(n)  str_replace(n,'\"','') %>% str_replace(.,'\"',''))
rownames(Data) <- otus
Tax[,1]<-otus

### Independent Variable
MAP <- read.csv('CP_map.txt',sep='\t',header = T)
X <- MAP$pH
names(X) <- MAP$X.SampleID
nms <- sapply(as.list(MAP$X.SampleID),toString)
nms <- sapply(as.list(nms),FUN=function(x) paste('X.',x,'.',sep=''))
Data <- Data[,nms]

### Phylogeny
tree <- read.tree('rep_set_aligned_gapfilt20_filt2.tre')
otus <- tree$tip.label
Data <- Data[tree$tip.label,]
Groups <- getGroups(tree)


# Preliminary Analysis ----------------------------------------------------

### Preliminary analysis - checking mean/variance scaling & P-value distributions
m <- colMeans(Data)
v <- apply(Data,MARGIN=2,var)
par(mfrow=c(1,1))
plot(m,v,log='xy')
lm(log(v)~log(m))

Dclr <- apply(Data,MARGIN=2,FUN=clr)
GLMs <- apply(Dclr,MARGIN=1,FUN=function(y,x) glm(y~x),x=X)
Pvals <- NULL
Pvals$PH <- sapply(GLMs,FUN=function(gg) summary(aov(gg))[[1]][1,'Pr(>F)'])
par(mfrow=c(2,2))
hist(Pvals$PH,main='PH Pvalue hist')

GLMs <- apply(Dclr,MARGIN=1,FUN=function(y,x) glm(y~x),x=MAP$Moisture)
Pvals$Moisture <- sapply(GLMs,FUN=function(gg) summary(aov(gg))[[1]][1,'Pr(>F)'])
hist(Pvals$Moisture,main='Moisture Pvalue hist')

GLMs <- apply(Dclr,MARGIN=1,FUN=function(y,x) glm(y~x),x=MAP$C)
Pvals$C <- sapply(GLMs,FUN=function(gg) summary(aov(gg))[[1]][1,'Pr(>F)'])
hist(Pvals$C,main='C Pvalue hist')


GLMs <- apply(Dclr,MARGIN=1,FUN=function(y,x) glm(y~x),x=MAP$N)
Pvals$N <- sapply(GLMs,FUN=function(gg) summary(aov(gg))[[1]][1,'Pr(>F)'])
hist(Pvals$N,main='N Pvalue hist')


# PhyloFactorization ------------------------------------------------------


system.time(PhyloFactor(Data,tree,X,choice = 'F',nfactors=1,ncores = 7,Grps=Groups,Pval.Cutoff = 0.1))
# 68.121 seconds for one factor. 

# load('Soil_microbiome_PF')
system.time(PhyloFactor(Data,tree,X,choice='var',Grps=Groups,nfactors=1,ncores=7))

PF.big.F <- PhyloFactor(Data,tree,X,choice = 'F',nfactors=1000,stop.early = T,ncores = 7,Grps=Groups)
# save(list=ls(),file='Soil_microbiome_PF')

PF.big <- PhyloFactor(Data,tree,X,stop.early = T,ncores = 7,Grps=Groups)
### NOTE: PF.big is JUST pH - for multiple regression used in Washburne et al. (2016) study, see PF.Mult below
# save(list=ls(),file='Documents/PhyloFactor/Soil_microbiome_PF_with_Tax')

ind <- which(!is.na(MAP$C))
log.C <- log(MAP$C[ind])
PF.C <- PhyloFactor(Data[,ind],tree,X=log.C,nfactors=5,ncores=7,Grps=Groups)



# load(file='Soil_microbiome_PF_with_Tax')
# ANALYSIS of PF.big ------------------------------------------------------

PF.big$nfactors #2091 factors
PF.big$factors[1:30,]

Cumulative.Ex.Var <- cumsum(PF.big$ExplainedVar)
# In total, regression of pH on the 2091 factors can only account for 12% of the variation in the data
# However, the first factor, which contains 1/6756 or 0.015% of the possible edges or 0.03% of the variables in the data,
# Captures 1.5% of the variation in the dataset. Although PF doesn't capture as much variation as PCA,
# it does identify edges that punch way above their weight in the dataset. 

phylofactor.visualize(PF.big,dimension=2,colorbar.name='pH')

smry <- phylofactor.summary(PF.big,taxonomy =Tax,factor=1)
td <- pf.tidy(smry)

# 206 monophyletic OTUs split from the rest
smry <- phylofactor.summary(PF.big,taxonomy =Tax,factor=1)
td <- pf.tidy(smry)
td$`group1, Monophyletic` ##Class Acidobacteria and class DA052 go down with pH compared to some other Acidobacteria and the rest
plot(sort(X),td$`Predicted ratio of group1/group2`[order(X)],type='l',lwd=4)
points(X,td$`Observed Ratio of group1/group2 geometric means`)

#193 monophyletic OTUs split from the rest
smry <- phylofactor.summary(PF.big,taxonomy =Tax,factor=2)
td <- pf.tidy(smry)
td$`group1, Monophyletic` ##Class Acidobacteria-6, [Cholarcidobacteria], an undescribed class of Acidobacteria, and c_S035 go up with pH
plot(sort(X),td$`Predicted ratio of group1/group2`[order(X)],type='l',lwd=4)
points(X,td$`Observed Ratio of group1/group2 geometric means`)

#31 monophyletic OTUs
smry <- phylofactor.summary(PF.big,taxonomy =Tax,factor=3)
td <- pf.tidy(smry)
td$`group2, Paraphyletic`
td$`group1, Monophyletic` ##an uncharacterized familiy of Actinomycetales and Actinomycetales f_Thermomonosporaceae is split from other Actinomycetales and the rest (including split from another uncharacterized family)
plot(sort(X),td$`Predicted ratio of group1/group2`[order(X)],type='l',lwd=4)
points(X,td$`Observed Ratio of group1/group2 geometric means`)
### These Actinomycetales go down with pH

#115 monophyletic OTUs
smry <- phylofactor.summary(PF.big,taxonomy =Tax,factor=4)
td <- pf.tidy(smry)
td$`group2, Paraphyletic`
td$`group1, Monophyletic` ##Two classes of Acidobacteria - Solibacteres and TM1 - go down with pH
plot(sort(X),td$`Predicted ratio of group1/group2`[order(X)],type='l',lwd=4)
points(X,td$`Observed Ratio of group1/group2 geometric means`)

#Factor 11 = 374 paraphyletic OTUs
smry <- phylofactor.summary(PF.big,taxonomy =Tax,factor=11)
td <- pf.tidy(smry)
td$`group2, Paraphyletic`
td$`group1, Paraphyletic` #A bunch of Bacteroidetes - c_Bacteroidia, c_Cytophagia, c_Flavobacteria, o_[Saprospirales], etc. - go up with pH.
plot(sort(X),td$`Predicted ratio of group1/group2`[order(X)],type='l',lwd=4,ylim=c(0.1,3))
points(X,td$`Observed Ratio of group1/group2 geometric means`)


# ANALYSIS of PF.C --------------------------------------------------------

PF.C$factors
phylofactor.visualize(PF.C,dimension=2,colorbar.name='C')
cumsum(PF.C$ExplainedVar)
## Carbon explains a very small fraction of the variance in the dataset.

FACTOR=1
smry <- phylofactor.summary(PF.C,taxonomy =Tax,factor=FACTOR)
td <- pf.tidy(smry)
td$`group1, Monophyletic` #The same as factor 1 for pH
plot(sort(log.C),td$`Predicted ratio of group1/group2`[order(log.C)],type='l',lwd=4,log='x')
points(log.C,td$`Observed Ratio of group1/group2 geometric means`)


FACTOR=2
smry <- phylofactor.summary(PF.C,taxonomy =Tax,factor=FACTOR)
td <- pf.tidy(smry)
td$`group1, Monophyletic` #The same as factor 1 for pH
plot(sort(log.C),td$`Predicted ratio of group1/group2`[order(log.C)],type='l',lwd=4,log='x')
points(log.C,td$`Observed Ratio of group1/group2 geometric means`)

FACTOR=3
smry <- phylofactor.summary(PF.C,taxonomy =Tax,factor=FACTOR)
td <- pf.tidy(smry)
td$`group1, Monophyletic` #The same as factor 1 for pH
plot(sort(log.C),td$`Predicted ratio of group1/group2`[order(log.C)],type='l',lwd=4,log='x')
points(log.C,td$`Observed Ratio of group1/group2 geometric means`)


# ANALYSIS of PF.N -------------------------------------------------------

N <- MAP$N
ind <- which(!is.na(N))
N <- N[ind]
par(mfrow=c(1,2))
hist(N,main='Distribution of N') #Note - the distribution of N is very skewed - this could complicate regression
hist(log(N),main='Distribution of log(N)') #that's a bit better - we'll use this

PF.N <- PhyloFactor(Data[,ind],tree,X=log(N),nfactors = 10,ncores=3)

save(list=ls(),file='Soil_PhyloFactorization')


PF.N$factors

FACTOR = 10
smry <- phylofactor.summary(taxonomy=Tax,factor=FACTOR,PF=PF.N)
TD <-  pf.tidy(smry)
TD$`group1`
par(mfrow=c(1,1))
plot(sort(log(N)),TD$`Observed Ratio of group1/group2 geometric means`[order(N)])
lines(sort(log(N)),TD$`Predicted ratio of group1/group2`[order(N)],lwd=2,col='red')
#Factor=3 has 6 OTUs of the genus Nitososphaera, with reliably low abundances at high N, and varied abundances at low N.
# Factor 5 has 12 members of Acidobacteria order Ellin6513
# Factor 6 has 87 members all of the Verrucoomicrobia class Spartobacteria, almost all of the family Chthoniobacteraceaae of which 
# Candidatus Xiphinematobacter is a member (known endosymbiont of plant-parasitic nematodes). Some were genus Chtoniobacter and others g__DA101.
# Factor 7 contains 8 members of the genus Kaistobacter f__Sphingomonodaceae
# Factor 8 contains Gemmatimonadetes (one unclassified order, and then o__C114 and o__N1423WL)
# Factor 9 contains 42 members of the Verrucomicrobia class Pedosphaerae
# Finally, factor 10 contains 100 microbes all of the class Thermoleophilia - some of the family Gaiellaceae and some of the order Solirubrobacterales
# These have the common pattern of reliably low/moderate abundances at high N, and variable to high abundances at low N


# ANALYSIS of multiple regression pH, C and N ----------------------------
indC <- which(!is.na(MAP$C))
indN <- which(!is.na(MAP$N))

ind <- intersect(indC,indN)

pH <- MAP$pH[ind]
C <- log(MAP$C)[ind]
N <- log(MAP$N)[ind]

X <- data.frame(pH,C,N)
frmla = Data~pH+C+N

PF.Mult <- PhyloFactor(Data[,ind],tree,X,frmla=frmla,nfactors=6,ncores=3)

ss <- lapply(PF.Mult$glms,FUN=function(gg) summary(aov(gg)))
library(relaimpo)
LMG <- lapply(PF.Mult$glms,FUN= function(gg) calc.relimp(gg,type='lmg'))
pH.imp <- sapply(LMG,FUN=function(rr) rr@lmg['pH']/sum(rr@lmg))
# Using 'lmg', pH accounts for 92.87%, 89.78%, and 92.94% of the first three factors

png('C:/Users/Big Alculus/Documents/Boulder/PhyloStats/Paper Figures/Soil_OV_pH.png',width=600,height=400)
phylofactor.visualize(PF.Mult,X=pH,dimension=3,colorbar.name = 'pH',colorbar.ticks=5)
dev.off()
png('C:/Users/Big Alculus/Documents/Boulder/PhyloStats/Paper Figures/Soil_OV_C.png',width=600,height=400)
phylofactor.visualize(PF.Mult,X=round(100*C)/100,,dimension=3,colorbar.name = 'C',colorbar.ticks = 4)
dev.off()
png('C:/Users/Big Alculus/Documents/Boulder/PhyloStats/Paper Figures/Soil_OV_N.png',width=600,height=400)
phylofactor.visualize(PF.Mult,X=round(100*N)/100,dimension=3,colorbar.name = 'N',colorbar.ticks=4)
dev.off()



### Now, let's visualize the PF.Mult vs. pH in ILR space AND let's bin the taxa and sort those bins
nfacts=4
binned_data_observed <- binProjection(PF.Mult,nfacts,Tax)
binned_data_predicted <- binProjection(PF.Mult,nfacts,Tax,prediction=T)
library(plotrix)

binned_data_observed$otus

phylofactor.visualize(PF.Mult,X=pH,dimension=3,colorbar.name = 'pH',colorbar.ticks=5)

LEGEND <- rev(names(binned_data_observed$otus))
LEGEND <- c('Acidophobic Acidobacteria','Solibacteres and TM1','Actinomycetales','Acidocaberiia and DA052','Remainder Bin')
cols <- rainbow(nfacts+1)
ix <- order(pH)
tiff('C:/Users/Big Alculus/Documents/Boulder/PhyloStats/Paper Figures/Soil_Bin_Plot2.png',width=1000,height=1000)
par(mfrow=c(2,1))
stackpoly(t(binned_data_observed$Data[,ix]),stack=T,xlab='pH',ylab='Relative Abundance',main='Observed BPUs',col=cols,xaxlab = pH[ix],axis4=F,cex.main=2,cex.lab=2)
stackpoly(t(binned_data_predicted$Data[,ix]),stack=T,xlab='pH',ylab='Relative Abundance',main='Predicted BPUs',col=cols,xaxlab = pH[ix],axis4=F,cex.main=2,cex.lab=2)
legend(170,0.8,legend=LEGEND,fill=rev(cols),cex=2)
dev.off()
# save(list=ls(),file = 'Soil_PF_R_Workspace')


# load('Soil_PF_R_Workspace')



# Summary of taxa split by each factor ------------------------------------------

FACTOR=1

PF.Mult$factors
TD <- phylofactor.summary(taxonomy=Tax,factor=FACTOR,PF=PF.Mult) %>% pf.tidy
TD$`group1, Monophyletic`

plot(pH,TD$`Observed Ratio of group1/group2 geometric means`)
#Group 1, which contains a monophyletic group of 206 OTUs, 
#Contains the c__Acidobacteriia and c__DA052, split from the rest.d
# This group goes down with pH

FACTOR=2

TD <- phylofactor.summary(taxonomy=Tax,factor=FACTOR,PF=PF.Mult) %>% pf.tidy
TD$`group1, Monophyletic`
smry <- phylofactor.summary(taxonomy=Tax,factor=FACTOR,PF=PF.Mult,simplify.TaxaIDs = F)
smry$group1$IDs

plot(pH,TD$`Observed Ratio of group1/group2 geometric means`)
#Group 1 in factor=2  contains 31 monophyletic OTUs, 
# all in the o__Actinomycetales - some in the family Thermomonosporaceae
# and all the rest being unclassified at the family level.
# These Actinomycetales go down with pH. ds


FACTOR=3

TD <- phylofactor.summary(taxonomy=Tax,factor=FACTOR,PF=PF.Mult) %>% pf.tidy
TD$`group1, Monophyletic`
plot(pH,TD$`Observed Ratio of group1/group2 geometric means`)
# The third factor splits more Acidobacteria - c__Solibacteres and c__TM1 - from the remainder.
# These Acidobacteria decrease in abundance with pH.



FACTOR=4
smry <- phylofactor.summary(taxonomy=Tax,factor=FACTOR,PF=PF.Mult)
TD <- phylofactor.summary(taxonomy=Tax,factor=FACTOR,PF=PF.Mult) %>% pf.tidy
TD$`group1, Monophyletic`
TD$`group2, Paraphyletic`
plot(pH,TD$`Observed Ratio of group1/group2 geometric means`)

## INCREASE with pH


FACTOR=5

TD <- phylofactor.summary(taxonomy=Tax,factor=FACTOR,PF=PF.Mult) %>% pf.tidy
TD$`group1, Monophyletic`
plot(pH,TD$`Observed Ratio of group1/group2 geometric means`)


FACTOR=6

TD <- phylofactor.summary(taxonomy=Tax,factor=FACTOR,PF=PF.Mult) %>% pf.tidy
TD$`group1, Monophyletic`
plot(pH,TD$`Observed Ratio of group1/group2 geometric means`)



############# Visualizing Tree #######################

atms <- bins(PF.Mult$basis[,1:3])
otus <- lapply(atms,FUN = function(x,names) names[x],names=tree$tip.label)
taxa <- lapply(otus,FUN=OTUtoTaxa,Taxonomy=Tax,common.name=F,uniques=T)

legendNames <- binTaxa(PF.Mult,3,Tax)
legendNames[[1]] <- 'Remainder'

tiff(filename = 'C:/Users/Big Alculus/Documents/Boulder/PhyloStats/Paper Figures/Soil_Phylogeny.tiff',width=2300,height=1500)
L = ColorTaxa(tree,Tax,type='unrooted',outputlegend = T,show.tip.label=F,scramble = T,use.edge.length=F)
lms <- par('usr')
legend(lms[1],lms[4],legend=L$Taxa,fill=L$colors,cex=2)
dev.off()

L = binPhyloPlot(PF.Mult,4,legend.return=T,edge.width=1,use.edge.length=F,type='unrooted')
lms=par('usr')
Legend <- lapply(L$bins,FUN = function(otus,tax,nms) OTUtoTaxa(nms[otus],tax,common.name=F,uniques = T),tax=Tax,nms=rownames(PF.Mult$Data))
legend(lms[1],lms[4],legend=Legend,fill=L$Colors)
