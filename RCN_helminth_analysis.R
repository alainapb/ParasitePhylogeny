## RUN THIS ONLY IF YOU ARE IN THE WORKING DIRECTORY
library(ape)
library(plyr)
library(tidyverse)
library(ouch)
library(dplyr)
library(phytools)
library(geiger)
library(caper)
library(picante)

## RCN_code_May2019 for small phylacine tree ##PLUS pruning out the 4 misplaced nematodes

## Read in GMPD data
GMPDmain<-read.csv("GMPD_main.csv")
GMPDparasite<-read.csv("GMPD_parasite_traits.csv")

## Read in the helminth tree
tree<-read.nexus("hel2_44burn_copy.trees")
## load the binomial names for the coded tree names
nom <- read.csv("name_key.csv")
## ensure that every tip has a binomial name
which(!(tree$tip.label %in% nom[,2])) ## ZTGLUN does not have a binomial name
## drop it from the tree
tree <- drop.tip(tree, "ZTGLUN")
tree <- drop.tip(tree, "HABMUS") #nematodes in the wrong area of the tree
tree <- drop.tip(tree, "ASCLUM") #nematodes in the wrong area of the tree
tree <- drop.tip(tree, "ASCSUU") #nematodes in the wrong area of the tree
tree <- drop.tip(tree, "PELSTR") #nematodes in the wrong area of the tree

write.tree(tree, file="~/Desktop/helminthtree.trees" )
write.tree(tree, file="~/Desktop/helminthtree.nex" )


## Load the generalism scores for every parasite
load("get_nriABUND.Rda") ## nri


########################## ANALYSES OF PHYLOGENETIC SIGNAL ##############################

#### GENERALISM AS MEAN PAIRWISE PHYLOGENETIC DISTANCE ####

## prune the tree to only include species with generalism scores
sapply(tree$tip.label, function(n) which(nom[,2]==n)) %>% unlist -> w
gen <- sapply(as.character(nom[w,1]), function(n) nri[rownames(nri)==n,"mpd.obs"]) %>% unlist
## turn all NA gen scores to zeros
gen[is.na(gen)] <- 0
## keep only the tips that have generalism scores
sapply(names(gen), function(n) as.character(nom[nom[,1]==n,2])) %>% unlist %>% unname -> keepers
mpd_tree <- drop.tip(tree, tree$tip.label[which(!(tree$tip.label%in%keepers))])
## relabel the tips with the binomial name
mpd_tree$tip.label <- names(gen)

## regular tree
contMap(mpd_tree, gen, fsize= .3, lwd=.75, outline=FALSE, sig = .5, leg.text="x")
add.scale.bar(length =.2)

## fan tree
contMap(mpd_tree, gen, fsize= .3, lwd=1.2, outline=FALSE, sig = .5, type = "fan")
add.scale.bar(length =.2)

## fit phylogenetic models to the data
fitContinuous(mpd_tree,gen, model="lambda") ##Best model: lambda=.5453; aicc=1861.819
fitContinuous(mpd_tree,gen,model="BM") #aicc=1943.006
fitContinuous(mpd_tree,gen,model="OU") #alpha = 2.71828; aicc=1905.292
fitContinuous(mpd_tree,gen,model="white") #aicc= 1887.075529
## other models: "trend","kappa","delta","drift", "EB"

#### GENERALISM AS NUMBER OF TAXA INFECTED ####

## prune the tree to only include species with number of host scores
sapply(tree$tip.label, function(n) which(nom[,2]==n)) %>% unlist -> j
no_taxa <- sapply(as.character(nom[j,1]), function(n) nri[rownames(nri)==n,"ntaxa"]) %>% unlist
## turn all NA ntaxa scores to zeros
no_taxa[is.na(no_taxa)] <- 0
## keep only the tips that have number of host scores
sapply(names(no_taxa), function(n) as.character(nom[nom[,1]==n,2])) %>% unlist %>% unname -> keepers
ntaxa_tree <- drop.tip(tree, tree$tip.label[which(!(tree$tip.label%in%keepers))])
## relabel the tips with the binomial name
ntaxa_tree$tip.label <- names(no_taxa)


### View tree with full names
plot(ntaxa_tree, cex=0.3)
add.scale.bar(length =.2)

## regular tree
contMap(ntaxa_tree, no_taxa, fsize= .3, lwd=.75, outline=FALSE, sig = .5, leg.text="x")
add.scale.bar(length =.2)

## fan tree
numbertree<-contMap(ntaxa_tree, no_taxa, fsize= .3, lwd=1.2, outline=FALSE, sig = .5, type = "fan")
add.scale.bar(length =.2)

## fit phylogenetic models to the data
fitContinuous(ntaxa_tree, log(no_taxa), model="lambda") #lambda=.2228; aicc=526.6889
fitContinuous(ntaxa_tree, log(no_taxa), model="BM") #aicc=675.1637
fitContinuous(ntaxa_tree, log(no_taxa), model="OU") #alpha=2.718; aicc=626.429
fitContinuous(ntaxa_tree, log(no_taxa), model="white") #aicc=527.850


##### CALCULATE MPD AMONG THE PARASITES INFECTING EACH HOST ######

## Load GMPD data and subset down to the species in the parasite phylogeny
GMPDsmall<-GMPDmain[GMPDmain$ParasiteCorrectedName %in% mpd_tree$tip.label, ]
GMPDsmall$HostCorrectedName <-gsub(" ", "_", GMPDsmall$HostCorrectedName)

## load mammal tree from Phylacine data
mam_tree<-read.nexus("mammal_phylo_consensus.nex") #load in consensus mammal tree to use as a response variable
small_mam_tree<-read.nexus("Small_phylogeny_consensus.nex") #load in consensus small mammal tree
#this tree is supposed to be more accurate, but may have less data (seems to have 2 less species of mammals)

##prune small tree to only include mammals from the GMPD
setdiff(small_mam_tree$tip.label, GMPDsmall$HostCorrectedName)->bye2
small_mam_tree<-drop.tip(small_mam_tree, small_mam_tree$tip.label[which(small_mam_tree$tip.label%in%bye2)])
plot(small_mam_tree, cex=0.3)

##calculate the mpd of parasites and hosts!
dis_para<-cophenetic.phylo(mpd_tree) #interspecific distance matrix of parasites
dis_mam<-cophenetic.phylo(small_mam_tree) #intespecific distance matrix of mammals from small tree


##make the pipeline to get the absense presence of each parasite (rows) that infect a given host (column)
mamlist<-as.vector(small_mam_tree$tip.label) #list of mammals in tip label order
para_list<-as.vector(mpd_tree$tip.label)
df <- data.frame(matrix(ncol=length(mamlist),nrow=length(para_list)), row.names = para_list)
colnames(df)<-mamlist

for (i in 1:nrow(df))
  df[i,] <- (colnames(df) %in% unique(subset(GMPDsmall, ParasiteCorrectedName==rownames(df)[i])$HostCorrectedName)) %>% as.numeric

write.csv(df,file = "df_test.csv", row.names = TRUE, col.names = TRUE)

## Create MPD scores for every parasite using the Phylacine tree (rather than ToL)
para.mpd <- mpd(df, dis_mam)
para.mpd[is.na(para.mpd)] <- 0 ## specialists have MPD of 0

## compare previously calculated parasite MPD to this new metric
merge(data.frame(parasite=names(gen), ToL.MPD=unname(gen)),
      data.frame(parasite=rownames(df), Phylacine.MPD=para.mpd),
      by="parasite") -> compMPD
with(compMPD, plot(ToL.MPD, Phylacine.MPD))


##check phylogenetic signal with the mpd created using the phylacine tree
#MAKE para.mpd a named number
names(para.mpd)<-para_list
#view phylacine mpd on tree
phyla_mpd_tree<-mpd_tree
contMap(phyla_mpd_tree, para.mpd, fsize= .3, lwd=1.0, outline=TRUE, sig = .5, leg.text="x", type="fan", legend = FALSE)
add.scale.bar(length=.5)

##test phylogenetic signal--has similar results to the ToL mpd scores
fitContinuous(phyla_mpd_tree, para.mpd, model="lambda") #lambda=.5336, aicc=1993.01
fitContinuous(phyla_mpd_tree, para.mpd, model="BM") #aicc=2060.50
fitContinuous(phyla_mpd_tree, para.mpd, model="OU") #alpha=2.718; aicc=2025.90
fitContinuous(phyla_mpd_tree, para.mpd, model="white") #aicc=2018.32



########################## PHYLOGENETIC REGRESSIONS ##############################

### Load host tree and host data

## adding in host trait data from Phylacine
Trait_data <- read_csv("PHYLACINE_DATA/MegaPast2Future-PHYLACINE_1.2-975f93c/Data/Traits/Trait_data.csv")

#checking the name consistency in the two dataframes
setdiff(GMPDsmall$HostCorrectedName, Trait_data$Binomial.1.2)
#correcting for spaces in the binomial names
GMPDsmall$HostCorrectedName<-gsub(" ", "_", GMPDsmall$HostCorrectedName)
#checking again for consistency in names between the two dataframes
setdiff(GMPDsmall$HostCorrectedName, Trait_data$Binomial.1.2)
#adding synomy data from phylacine to use the Elton binomial (matches GMPD)
syn <- read_csv("PHYLACINE_DATA/MegaPast2Future-PHYLACINE_1.2-975f93c/Data/Taxonomy/Synonymy_table_valid_species_only.csv")
#taxonomy is split in to two columns, merge them together
syn$Eltonbinomial<-paste0(syn$EltonTraits.1.0.Genus, syn$EltonTraits.1.0.Species)
#change spaces in binomial name to "_"
syn$Eltonbinomial<-paste0(syn$EltonTraits.1.0.Genus,"_", syn$EltonTraits.1.0.Species)
#check differences between naming
setdiff(GMPDsmall$HostCorrectedName, syn$Eltonbinomial)
#change names to be consistent...
GMPDsmall$HostCorrectedName<-gsub("Cercopithecus_albogularis", "Cercopithecus_mitis", GMPDsmall$HostCorrectedName)
setdiff(GMPDsmall$HostCorrectedName, syn$Eltonbinomial)
##change 'Equus_burchelli' to 'Equus_??' in GMPDsmall
GMPDsmall$HostCorrectedName<-gsub("Equus_burchelli", "Equus_quagga", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Equus_quaggai", "Equus_quagga", GMPDsmall$HostCorrectedName)
##remove the ones without a name
GMPDsmall <- GMPDsmall[GMPDsmall$HostCorrectedName!="no_binomial_name",]

### Load parasite tree and parasite data
## Subset the parasite transmission mode data to the parasites in the tree
GMPDparasite<-GMPDparasite[GMPDparasite$CorrectName %in% mpd_tree$tip.label, ]

## Load parasite body size data from Benesh, Lafferty, and Kuris 2017
## Length and width measurments are in mm 
benesh_data_para<- read_csv("Benesh_et_al_2017/ecy1680-sup-0001-datas1/CLC_database_lifehistory.csv")
benesh_data_para$Host.species <-gsub(" ", "_", benesh_data_para$Host.species)

## remove parasites from Benesh et al. dataset that don't exist in the phylogeny
setdiff(benesh_data_para$Parasite.species, GMPDsmall$ParasiteCorrectedName)->nothanks
benesh_data_para[-which(benesh_data_para$Parasite.species %in% nothanks) ,] -> b_para_subset

## average the body size for rows that are the same species
subset(b_para_subset, Sex=="f") %>%
  group_by(Parasite.species) %>%
  summarise(mean_l=mean(Length, na.rm = TRUE)) ->data_fem_len

subset(b_para_subset, Sex=="f") %>%
  group_by(Parasite.species) %>%
  summarise(mean_w=mean(Width, na.rm = TRUE)) ->data_fem_wid

## combine length and width in a meaningful way to calculate 'body size' (L^2*W)?
b.size<-left_join(data_fem_len, data_fem_wid, by="Parasite.species")
b.size$ratio<-(b.size$mean_l*b.size$mean_w)

write.csv(b.size, "bodysizepara.csv")


##Run same predictor variables on 3 responses: 1) parasite mpd 2) parasite host #

## SUBSET DATA FOR PGLS ##
## Trim the tree to include only those parasite species for which
## there is body size data
setdiff(mpd_tree$tip.label, b.size$Parasite.species)->no_body_size
phy <- drop.tip(mpd_tree, mpd_tree$tip.label[which(mpd_tree$tip.label%in%no_body_size)])
## make sure body size is in the same order as the parasite phylogeny
b.size <- b.size[match(phy$tip.label, b.size$Parasite.species),]

## Data on number of taxa each parasite infects for parasites for which
## there is body size data
no_taxa1 <- as.data.frame(no_taxa)
no_taxa1$Para <- ntaxa_tree$tip.label
setdiff(no_taxa1$Para, b.size$Parasite.species)->no_body_size1
no_taxa1[-which(no_taxa1$Para %in% no_body_size1) ,] -> no_taxa_subset
## make sure this is in the same order as the parasite phylogeny
phy$tip.label==no_taxa_subset$Para

## Data on the MPD for each parasite for which there is body size data
gen1<-as.data.frame(gen)
gen1$Para<-mpd_tree$tip.label
setdiff(gen1$Para, b.size$Parasite.species)->no_body_size2
gen1[-which(gen1$Para %in% no_body_size2) ,]->gen_subset
## make sure this is in the same order as the parasite phylogeny
phy$tip.label==gen_subset$Para

## Data on the transmission mode for each parasite
setdiff(GMPDparasite$CorrectName, b.size$Parasite.species)->no_body_size3
GMPDparasite[-which(GMPDparasite$CorrectName %in% no_body_size3) ,] -> GMPDparasite_subset
## make sure this is in the same order as the parasite phylogeny
phy$tip.label==GMPDparasite_subset$CorrectName
## resort to get in the right order
GMPDparasite_subset <- GMPDparasite_subset[match(phy$tip.label, GMPDparasite_subset$CorrectName),]

## RUNNING PGLS ##
## Combine all data (body size as length, number of infected taxa, MPD, and transmission mode) into a single data.frame
dat<-data.frame(taxa=phy$tip.label,
                para.l=b.size$mean_l,
                para.w=b.size$mean_w,
                para.b=b.size$ratio,
                no.taxa=no_taxa_subset$no_taxa,
                para.mpd=gen_subset$gen,
                close=as.numeric(GMPDparasite_subset$close),
                nonclose=as.numeric(GMPDparasite_subset$nonclose),
                intermediate=as.numeric(GMPDparasite_subset$intermediate))

## seems like there might be a weak relationship between MPD and parasite body size
with(dat, plot(para.l, para.mpd))

## doesn't seem like there's likely to be anything here
dat %>% group_by(close) %>% summarise(mean.mpd = mean(para.mpd),
                                      sd.mpd=sd(para.mpd),
                                      mean.ntaxa = mean(no.taxa),
                                      sd.ntaxa=sd(no.taxa))
dat %>% group_by(nonclose) %>% summarise(mean.mpd = mean(para.mpd),
                                         sd.mpd=sd(para.mpd),
                                         mean.ntaxa = mean(no.taxa),
                                         sd.ntaxa=sd(no.taxa))
dat %>% group_by(intermediate) %>% summarise(mean.mpd = mean(para.mpd),
                                             sd.mpd=sd(para.mpd),
                                             mean.ntaxa = mean(no.taxa),
                                             sd.ntaxa=sd(no.taxa))


## constructs the data.frame of the taxa names and data in the order of the tree tips
cdat<-comparative.data(data=dat, phy=phy, names.col = "taxa")
print(cdat)

## model testing
## this works! But no aspect of body size has any affect on MPD in this regression
mod0 <- pgls(para.mpd~para.l, cdat, lambda = "ML")
anova(mod0)
summary(mod0)

mod00<-pgls(para.mpd~para.l, cdat, lambda = .55)
anova(mod00)
summary(mod00)

mod01 <- pgls(para.mpd~para.w, cdat, lambda = "ML")
anova(mod01)

mod001 <- pgls(para.mpd~para.w, cdat, lambda = .55)
anova(mod001)

mod02 <- pgls(para.mpd~para.b, cdat, lambda = "ML")
anova(mod02)

mod002 <- pgls(para.mpd~para.b, cdat, lambda = .55)
anova(mod002)

###PGLS on mpd
## this does not work
mod1<-pgls(para.mpd ~ close * nonclose * intermediate, cdat, lambda = 'ML')

## this does work, suggesting the issue is dealing with interactions 
mod1<-pgls(para.mpd ~ close + nonclose + intermediate + para.l , cdat, lambda = 'ML')
summary(mod1)
plot(mod1)

mod1_<-pgls(para.mpd ~ close + nonclose + intermediate + para.l , cdat, lambda = .55)
summary(mod1_)
plot(mod1_)

mod1.s<-pgls(para.mpd ~ close + nonclose + intermediate, cdat, lambda = .55)
summary(mod1.s)

#### PGLS on No. of host taxa

mod2<-pgls(no.taxa ~ close + nonclose + intermediate + para.l, cdat, lambda = 'ML')
summary(mod2)
plot(mod2)

mod2_<-pgls(no.taxa ~ close + nonclose + intermediate + para.l, cdat, lambda = .22)
summary(mod2_)
plot(mod2_)

mod2_s<-pgls(no.taxa ~ para.b, cdat, lambda = .22)
summary(mod2_s)

mod2.s<-pgls(no.taxa ~ close + nonclose + intermediate, cdat, lambda = .22)
summary(mod2.s)


##Generalism scores for mammals to plot on to mammal tree
##Use to identify clades of mammals with phylogenetically diverse parasites..?

## Create MPD scores for every host using the helminth tree
mam.mpd <- mpd(t(df), dis_para)
mam.mpd[is.na(mam.mpd)] <- 0 ## hosts infected by only one parasite have MPD of 0

#add names to mam.mpd
mamlist<-gsub(" ", "_", mamlist)
names(mam.mpd)<-mamlist

## regular tree
contMap(small_mam_tree, mam.mpd, fsize= .3, lwd=.75, outline=FALSE, sig = .5, leg.text="x")
add.scale.bar(length =50)

## fan tree
contMap(small_mam_tree, mam.mpd, fsize= .35, lwd=1.2, outline=TRUE, sig = .5, type = "fan", legend = FALSE)
add.scale.bar(length= 50)

## fit phylogenetic models to the data
fitContinuous(small_mam_tree, mam.mpd, model="lambda") #


fitContinuous(small_mam_tree, mam.mpd,model="BM") #
fitContinuous(small_mam_tree, mam.mpd,model="OU") #
fitContinuous(small_mam_tree, mam.mpd,model="white") #
## other models: "trend","kappa","delta","drift", "EB"


#####Detection Bias??#####

##test data bias, correlation between para.mpd and no. of citations
a<-as.data.frame(table(GMPDsmall$ParasiteCorrectedName))
setdiff(a$Var1, para_list)->kickout
a[-which(a$Var1 %in% kickout),]->a

freq<-a[match(phyla_mpd_tree$tip.label, a$Var1),]

cor(a$Freq, para.mpd) #r=-0.0398
ggplot(a, mapping = aes(a$Freq,para.mpd)) + geom_point() ##Echinococcus multilocularis has 429 observations in the dataset

##without E. multilocularis --no obvious bias
ggplot(a, mapping = aes(a$Freq,para.mpd)) + 
  geom_point() +
  xlim(0,150)+
  geom_smooth(method='lm', col='grey38')+
  labs(x='No. of Parasite Observations', y='Parasite MPD')+
  theme_classic()

##test mammal bias
b<-as.data.frame(table(GMPDsmall$HostCorrectedName))
setdiff(b$Var1, mamlist)->kickout2
b[-which(b$Var1 %in% kickout2),]->b

freq2<-b[match(small_mam_tree$tip.label, b$Var1),]

cor(b$Freq, mam.mpd) #r=.0592
ggplot(b, mapping = aes(b$Freq, mam.mpd)) + geom_point() #Vulpes vulpes has 927 observations

##without Vulpes vulpes
ggplot(b, mapping = aes(b$Freq, mam.mpd)) + 
  geom_point() +
  xlim(0,220)+
  geom_smooth(method='lm', col='grey38')+
  labs(x='No. of Mammal Observations', y='Mammal MPD')+
  theme_classic()

