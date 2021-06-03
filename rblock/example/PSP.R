library(igraph)
library(dplyr)
library(rblock)


runModel <- function(GN, LABELS, C=-1, ITS=10, TITLE="", MODEL=c(0), NCORES=0){

    RES <- list()

    if( C == -1 ){    
        C    = length(table(LABELS))
    } 
        
    anno = cbind(GN,LABELS)
    anno = as.data.frame(anno)

    print(PRINT=1)
    loadMetaData(DF=anno,MODELS=MODEL)

    runDCSBM(K=C,SLOT=0,ITS=ITS,NCORES=NCORES)
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

## load network
gg = read.graph("PPI_PSP.gml", format="gml")

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

#--choose which study to run on PSP network
runStudies <- vector(length=6)
runStudies[1] = 0
runStudies[2] = 0
runStudies[3] = 0
runStudies[4] = 0
runStudies[5] = 1
runStudies[6] = 0

#--container to store results
STUDY <- list()

if( runStudies[1] ){

## Study 1, AD, PD and HTN
title = "AD_PD_HTN"
STUDY[[1]] = runModel(GN=V(gg)$name,LABELS=label3,ITS=10,TITLE=title)
names(STUDY)[1] = title

}

if( runStudies[2] ){

## Study 2, AD and HTN
title = "AD_HTN"
STUDY[[2]] = runModel(GN=V(gg)$name,LABELS=label1,ITS=10,TITLE=title)
names(STUDY)[2] = title

}

if( runStudies[3] ){

## Study 3, AD and HTN
title = "PD_HTN"
STUDY[[3]] = runModel(GN=V(gg)$name,LABELS=label2,ITS=10,TITLE=title)
names(STUDY)[3] = title

}

if( runStudies[4] ){

## Study 4, AD and PD
title = "AD_PD"
STUDY[[4]] = runModel(GN=V(gg)$name,LABELS=label4,ITS=10,TITLE=title)
names(STUDY)[4] = title

}

if( runStudies[5] ){

## Study 5, AD and HTN
title = "AD_HTN_46"
STUDY[[5]] = runModel(GN=V(gg)$name,LABELS=label1,C=46,ITS=10,TITLE=title,NCORES=10)
names(STUDY)[5] = title

}


if( runStudies[6] ){

## Study 6, No diseases
title = "None_46"
STUDY[[6]] = runModel(GN=V(gg)$name,LABELS=non,C=46,ITS=10,TITLE=title,NCORES=10)
names(STUDY)[6] = title

}
    
saveRDS(STUDY, "SBMstudies.RDS",compress=TRUE)
