library(igraph)
library(dplyr)
library(rblock)

addNoise <- function( Y, MN=0, SD=0.05 ){
    return( Y+rnorm(length(Y),mean=MN, sd=SD) )
}

scale <- function(X){
    
    X = as.numeric(as.vector(X))

    xmin <- min(X,na.rm=T)
    xmax <- max(X,na.rm=T)

    return((X-min(X))/(max(X)-min(X)))
    
}

getCDF <- function(X){

    X = as.numeric(as.vector(X))
    
    mn = floor(min(X))
    mx = ceiling(max(X))
    #mx = floor(max(X))

    cdf = vector(length=length(X))

    z = seq(mn,mx,by=0.01)
    P = ecdf(X)
    p=P(z)

    for(i in 1:length(X)){
        cdf[i] = last(p[z <= X[i]])
    }

    return(cdf)

}

unitvec <- function(X){
    X = as.numeric(as.vector(X))
    return( X/(sqrt(sum(X^2))) )    
}

tt = read.delim("directed.csv",sep="\t",header=F)
tt = as.data.frame(tt)
print(PRINT=1)
loadGraph(DF=tt,directed=1)

#dir="testData"
#dir=""

#gg = read.graph(sprintf("%s/NRCMASC.gml",dir), format="gml")
gg = read.graph("NRCMASC.gml", format="gml")

ed = get.edgelist(gg,names=T)
ED = as.data.frame(ed)

loadGraph(DF=ED)

run=1

if( run ){
    oo    = cbind(V(gg)$name,V(gg)$SpecMod)
    test  = addNoise(as.numeric(oo[,2]),MN=0,SD=1)       
}

anno   = cbind(V(gg)$name, V(gg)$Disease, test)
anno   = as.data.frame(anno)
models = c(0,1)

print(PRINT=1)
loadMetaData(DF=anno,MODELS=models)
runDCSBM(K=10,SLOT=1)
getVertexKProbs()


run=0

if( run ){

    oo = cbind(V(gg)$name,V(gg)$SpecMod)
    x  = addNoise(as.numeric(oo[,2]),MN=0,SD=1)
    oo = cbind(oo,x)

    
}
