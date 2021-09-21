library(ape)
library(tidyverse)
library(phytools)
library(geiger)

## load tree
tree<-read.nexus("hel2_44burn_copy.trees")
## load the binomial names for the coded tree names
nom <- read.csv("name_key.csv")
# make nom characters not factors
nom$CorrectName<-as.character(nom$CorrectName)
nom$Code<-as.character(nom$Code)
## ensure that every tip has a binomial name
which(!(tree$tip.label %in% nom[,2])) ## ZTGLUN does not have a binomial name - code should actually be ZYGLUN; fix
tree$tip.label[which(!(tree$tip.label %in% nom[,2]))] <- "ZYGLUN"

## drop 5 nematodes and a acanth that were clearly placed in the wrong areas of the tree
tree <- drop.tip(tree, c("HABMUS","ASCLUM","ASCSUU","PELSTR","FILMAR","PROELE")) # 

## binomial names for only the species in the tree
nom[nom$Code %in% tree$tip.label,]->a

##nom in same order as tree tip labels
b<- a[match(tree$tip.label, a$Code),]

#replace code names with binomial names
tree$tip.label<-b$CorrectName

#rename the global tree for GMPD data
global_tree<-tree

## Read in GMPD data
GMPDmain<-read.csv("GMPD_main.csv")
GMPDparasite<-read.csv("GMPD_parasite_traits.csv")

#host-parasite association data
nearctic_assoc_data<-read.csv("Host_Range_Nearctic_Mammals_list.csv", header=T)

#change to characters instead of factors
nearctic_assoc_data$Host<-as.character(nearctic_assoc_data$Host)
nearctic_assoc_data$Parasite<-as.character(nearctic_assoc_data$Parasite)

## Load GMPD data and subset down to the species in the parasite phylogeny
GMPDsmall<-GMPDmain[GMPDmain$ParasiteCorrectedName %in% global_tree$tip.label, ]
GMPDsmall[which(GMPDsmall$HostReportedName=="Macaca hecki / M. tonkeana hybrid"),"HostCorrectedName"] <- "Macaca hecki" ## based on details from that study
GMPDsmall[which(GMPDsmall$HostReportedName=="Hapalemur sp."),"HostCorrectedName"] <- "Hapalemur griseus" ## based on collection location, this is the only species found there
GMPDsmall[which(rownames(GMPDsmall)=="11772"),"HostCorrectedName"] <- "Papio anubis" ## based on collection location
GMPDsmall[which(rownames(GMPDsmall)=="11773"),"HostCorrectedName"] <- "Papio ursinus" ## based on collection location
GMPDsmall[which(rownames(GMPDsmall)=="11862"),"HostCorrectedName"] <- "Papio anubis" ## based on collection location
GMPDsmall[which(rownames(GMPDsmall)=="11863"),"HostCorrectedName"] <- "Papio anubis" ## based on collection location
GMPDsmall[which(rownames(GMPDsmall)=="11864"),"HostCorrectedName"] <- "Papio anubis" ## based on collection location
GMPDsmall[which(rownames(GMPDsmall)=="11865"),"HostCorrectedName"] <- "Papio ursinus" ## based on collection location
GMPDsmall[which(rownames(GMPDsmall)=="11866"),"HostCorrectedName"] <- "Papio ursinus" ## based on collection location
GMPDsmall[which(rownames(GMPDsmall)=="11903"),"HostCorrectedName"] <- "Papio anubis" ## based on collection location

##load complete mammal tree from PHYLACINE
large_mam_tree<-read.nexus("mammal_phylo_consensus.nex")
large_mam_tree$tip.label<-gsub("_", " ", large_mam_tree$tip.label)

#### MAKE A MAMMAL TREE FOR HOSTS IN THE GMPD DATASET ####
## Make GMPD names match Elton binomial names found in the tree
GMPDsmall$HostCorrectedName<-gsub("Cercopithecus albogularis", "Cercopithecus mitis", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Equus burchelli", "Equus quagga", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Equus quaggai", "Equus quagga", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Otaria flavescens", "Otaria byronia", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Cercopithecus lhoesti","Allochrocebus lhoesti", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Equus caballus","Equus ferus", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Lama glama","Lama guanicoe", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Bos frontalis","Bos gaurus", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Equus asinus", "Equus hemionus", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Equus asinus", "Equus hemionus", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Monachus schauinslandi", "Neomonachus schauinslandi", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Neotragus moschatus", "Neotragus pygmaeus", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Ovis aries", "Ovis orientalis", GMPDsmall$HostCorrectedName)
GMPDsmall$HostCorrectedName<-gsub("Taurotragus oryx", "Tragelaphus oryx", GMPDsmall$HostCorrectedName)
##remove records from the GMPD (all from primates) where hosts were identified only to genus (super-weird that this happened!)
GMPDsmall <- GMPDsmall[GMPDsmall$HostCorrectedName!="no binomial name",]

##prune large tree to only include mammals from the GMPD
GMPD_mammal_tree<-drop.tip(large_mam_tree, setdiff(large_mam_tree$tip.label, GMPDsmall$HostCorrectedName))

#### MAKE A MAMMAL TREE FOR HOSTS IN THE NEARCTIC DATASET ####
## Modify the names in the Nearctic dataset to match the names on the large mammal tree
gsub("Mustela vison", "Neovison vison", nearctic_assoc_data$Host) -> nearctic_assoc_data$Host
## rename Alces alces (Eurasian elk) to Alces americanus (American moose) - these are sister species (or even subspecies) so any phylogenetic conclusions are unaltered in an analysis containing only one of the two species
gsub("Alces americanus", "Alces alces", nearctic_assoc_data$Host) -> nearctic_assoc_data$Host
## rename Arctic fox in the data
gsub("Alopex lagopus", "Vulpes lagopus", nearctic_assoc_data$Host) -> nearctic_assoc_data$Host
## rename American porcupine in the data
gsub("Erethizon dorsatus", "Erethizon dorsatum", nearctic_assoc_data$Host) -> nearctic_assoc_data$Host

## Remove Ascaris lumbricoides from the Nearctic dataset because it is a human parasite
nearctic_assoc_data <- nearctic_assoc_data[-which(nearctic_assoc_data$Parasite=="Ascaris lumbricoides"),]

##prune large tree to only include mammals from the Nearctic dataset
nearctic_mammal_tree<-drop.tip(large_mam_tree, setdiff(large_mam_tree$tip.label, nearctic_assoc_data$Host)) 
##subset global tree to nearctic data
nearctic_tree<-drop.tip(global_tree, setdiff(global_tree$tip.label, nearctic_assoc_data$Parasite))


#parasite abundance/occurence data
nearctic_abund_data<-read.csv("Host_Range_Nearctic_Mammals.csv")

#change to characters instead of factors
nearctic_abund_data$Parasite<-as.character(nearctic_abund_data$Var1)


## Need to control for sampling effort in future analyses





## save the relevant databases and trees
saveRDS(nearctic_mammal_tree, "Nearctic_mammal_tree.RDS")
saveRDS(nearctic_tree, "Nearctic_parasite_tree.RDS")
saveRDS(nearctic_assoc_data, "Nearctic_data.RDS")

saveRDS(GMPD_mammal_tree, "GMPD_mammal_tree.RDS")
saveRDS(global_tree, "GMPD_parasite_tree.RDS")
saveRDS(GMPDsmall, "GMPD_data.RDS")

## Idiot check: are all parasites and hosts found in both the trees and the association data?
Nearctic_mammal_tree <- readRDS("Nearctic_mammal_tree.RDS")
Nearctic_parasite_tree <- readRDS("Nearctic_parasite_tree.RDS")
Nearctic_data <- readRDS("Nearctic_data.RDS")

GMPD_mammal_tree <- readRDS("GMPD_mammal_tree.RDS")
GMPD_parasite_tree <- readRDS("GMPD_parasite_tree.RDS")
GMPD_data <- readRDS("GMPD_data.RDS")

setdiff(Nearctic_mammal_tree$tip.label, Nearctic_data$Host)
setdiff(Nearctic_data$Host, Nearctic_mammal_tree$tip.label)
setdiff(Nearctic_parasite_tree$tip.label, Nearctic_data$Parasite)
setdiff(Nearctic_data$Parasite, Nearctic_parasite_tree$tip.label)

setdiff(GMPD_mammal_tree$tip.label, GMPD_data$HostCorrectedName)
setdiff(GMPD_data$HostCorrectedName, GMPD_mammal_tree$tip.label)
setdiff(GMPD_parasite_tree$tip.label, GMPD_data$ParasiteCorrectedName)
setdiff(GMPD_data$ParasiteCorrectedName, GMPD_parasite_tree$tip.label)


