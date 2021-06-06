## clean
rm(list=ls())

library(igraph)
library(dplyr)
library(rblock)


runModel <- function(GN, LABELS, C=-1, ITS=10, TITLE="", MODEL=c(0), SEED=0, NCORES=0){

    RES <- list()

    if( C == -1 ){    
        C    = length(table(LABELS))
    } 
        
    anno = cbind(GN,LABELS)
    anno = as.data.frame(anno)

    print(PRINT=1)
    loadMetaData(DF=anno,MODELS=MODEL)

    runDCSBM(K=C,SLOT=0,ITS=ITS,SEED=SEED,NCORES=NCORES)
    res = getVertexKProbs()

    df   = as.data.frame(res$df)
    write.table(df,sprintf("study_df_%s.csv",TITLE),sep="\t",row.names=F,col.names=T,quote=F)

    np   = as.data.frame(res$np)
    colnames(np) = c("id",seq(0,(ncol(np)-2),1))
    write.table(np,sprintf("study_np_%s.csv",TITLE),sep="\t",row.names=F,col.names=T,quote=F)

    Vanno = as.data.frame(res$anno)
    write.table(Vanno,sprintf("study_Vanno_%s.csv",TITLE),sep="\t",row.names=F,col.names=T,quote=F)

    colnames(anno) = c("name","label")
    write.table(anno,sprintf("study_anno_%s.csv",TITLE),sep="\t",row.names=F,col.names=T,quote=F)

    
    RES$DF    = df
    RES$NP    = np
    RES$VANNO = Vanno
    RES$ANNO  = anno
    
    return(RES)
    
}

args <- commandArgs(TRUE);

SEED   <- as.numeric(args[1])
ITS    <- as.numeric(args[2])
CORES  <- as.numeric(args[3])
NTWRK  <- as.character(args[4])
K      <- as.numeric(args[5])

cat("SEED:  ", SEED, "\n")
cat("ITS:   ", ITS, "\n")
cat("CORES: ", CORES, "\n")
cat("NTWRK: ", NTWRK, "\n")
cat("K:     ", K, "\n")

#---Data directory
DIR  <- "datasets/"

#---Load networks
graphs  <- list.files(DIR,"*.gml")
Ngraphs <- unlist(strsplit(graphs,".gml"))

indx <- which(NTWRK==Ngraphs)

if( length(indx) != 0 ){

## load network
gg = read.graph(sprintf("%s/%s",DIR,graphs[indx]), format="gml")

## build edge list
ed = get.edgelist(gg,names=T)
ED = as.data.frame(ed)

loadGraph(DF=ED)

## build annotation set
dis = V(gg)$TopOntoOVG
non = rep("",length(dis))
ad  = grepl("AD",dis)
htn = grepl("HTN",dis)
pd  = grepl("PD",dis)

adNothtn       = ad  & !htn
htnNotad       = htn & !ad
adANDhtn       = ad  & htn

pdNothtn       = pd  & !htn
htnNotpd       = htn & !pd
pdANDhtn       = pd  & htn

adNotpd        = ad  & !pd
pdNotad        = pd  & !ad
adANDpd        = ad  &  pd

adNotpdNothtn =  ad & !pd & !htn
adANDpdNothtn =  ad &  pd & !htn
adNotpdANDhtn =  ad & !pd &  htn
adANDpdANDhtn =  ad &  pd &  htn
pdNotadNothtn = !ad &  pd & !htn
pdNotadANDhtn = !ad &  pd &  htn
htnNotadNotpd = !ad & !pd &  htn

label1 = ifelse(adNothtn,"AD","other")
label1 = ifelse(htnNotad,"HTN",label1)
label1 = ifelse(adANDhtn,"AD&HTN",label1)

label2 = ifelse(pdNothtn,"PD","other")
label2 = ifelse(htnNotpd,"HTN",label2)
label2 = ifelse(pdANDhtn,"PD&HTN",label2)

label3 = ifelse(adNotpdNothtn,"AD","other")
label3 = ifelse(adANDpdNothtn,"AD&PD",label3)
label3 = ifelse(adNotpdANDhtn,"AD&HTN",label3)
label3 = ifelse(adANDpdANDhtn,"AD&PD&HTN",label3)
label3 = ifelse(pdNotadNothtn,"PD",label3)
label3 = ifelse(pdNotadANDhtn,"PD&HTN",label3)
label3 = ifelse(htnNotadNotpd,"HTN",label3)

label4 = ifelse(adNotpd,"AD","other")
label4 = ifelse(pdNotad,"PD",label4)
label4 = ifelse(adANDpd,"AD&PD",label4)


#----

STUDY <- list()

## Study 1, AD and HTN
#title = sprintf("AD_HTN_%g",K)
#cat("run model: ", title,"\n")
#STUDY[[1]] = runModel(GN=V(gg)$name,LABELS=label1,C=K,ITS=ITS,TITLE=title,SEED=SEED,NCORES=CORES)
#names(STUDY)[1] = title

## Study 2, No diseases
title = "None_46"
STUDY[[1]] = runModel(GN=V(gg)$name,LABELS=non,C=K,ITS=ITS,TITLE=title,SEED=SEED,NCORES=CORES)
names(STUDY)[1] = title

saveRDS(STUDY, "SBMstudies.RDS",compress=TRUE)

} else {
  cat("Can't find graph!")
}