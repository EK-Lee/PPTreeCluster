#' @title Construct The Projection Pursuit Cluster Tree
#' @description  Find tree structure using various projection pursuit indices for
#' clustering in each split.
#' @usage PPTreeCluster(origdata,
#'                      PPmethod="Holes",
#'                      gamma=0.5,
#'                      method="kNN",
#'                      eps=0, 
#'                      prob=0.2,
#'                      weight=FALSE,
#'                      nk=0L,
#'                      t=1,
#'                      m=0L,
#'                      probT=0.05, 
#'                      R=0.1, 
#'                      max_depth=5,
#'                      minNodeN=5,
#'                      energy=0,
#'                      cooling=0.999,
#'                      TOL=1/1000000,
#'                      maxiter=50000,
#'                      seed=NULL,
#'                      notenough.plot = FALSE,
#'                      relCUTOFF=0.1,
#'                      cutB=0.2, 
#'                      minNrate=0.1,
#'                      ...)
#' @param origdata data matrix without class information
#' @param PPmethod projection pursuit method for clustering; "Holes", "Skew", 
#'                 "Comb", "Lpp", "NatrualHermite", and "P"
#' @param gamma combine ratio
#' @param method "kNN" determines neighbors based on count, while "dist" 
#'               determines based on distance.
#' @param eps distance criterion at the "dist" method
#' @param prob ratio for determining the neighborhood for kNN method        
#' @param weight weight flag in LDA, PDA and Lr index
#' @param nk count criterion at the "kNN" method
#' @param t weight parameter
#' @param m order of the index. low order captures the overall structure 
#'          of data, high order finds detailed patterns of data
#' @param probT the ratio to trim data
#' @param R criterion for determining neighbors
#' @param max_depth maximum depth of the projection pursuit cluster tree
#' @param minNodeN minimum number of observations in the final node
#' @param energy parameter for optimization
#' @param cooling parameter for optimization
#' @param TOL parameter for optimization
#' @param maxiter parameter for optimization
#' @param seed seed for random number
#' @param relCUTOFF parameter for density estimation
#' @param notenough.plot flag to draw plot
#' @param cutB parameter to define cluster
#' @param minNrate parameter to define cluster
#' @param ... Arguments to be passed to methods.
#' @return PPtreecluster object
#' @export
#' @examples
#' data(sampledata)
#' Cluster.result <- PPTreeCluster(sampledata[,-1])
#' Cluster.result

PPTreeCluster<-function(origdata,
                        PPmethod="Holes",
                        gamma=0.5,
                        method="kNN",eps=0, 
                        prob=0.2,weight=FALSE,nk=0L,t=1,m=0L,
                        probT=0.05, R=0.1, 
                        max_depth=5,
                        minNodeN=5,
                        energy=0,
                        cooling=0.999,
                        TOL=1/1000000,
                        maxiter=50000,
                        seed=NULL,
                        notenough.plot = FALSE,
                        relCUTOFF=0.1,
                        cutB=0.2, 
                        minNrate=0.1,...){
  origdata = as.matrix(origdata)
  if(!is.null(seed))
      set.seed(seed)
  if(is.null(minNodeN)){
    minNodeN = max(floor(nrow(origdata)*minNrate),10)
  }  else{
    minNodeN = max(floor(nrow(origdata)*minNrate),minNodeN)
  }

  ### Find b
  FindB <- function(proj.data,cutB=0.2){
    #plot(density(proj.data))
    hp0 <- mdh(proj.data)
    hp_new <-hp <- hp0
    for(i in 1:5){
      hp <- mdh(proj.data, v0 = hp$v, 
                bandwidth = hp$params$h*0.9, 
                alphamin = hp$params$alpha)
      if(abs(hp$rel.dep)>cutB & abs(hp$rel.dep)>abs(hp_new$rel.dep)) 
        hp_new <- hp
    }
    if(abs(hp$rel.dep)<cutB)
        hp_new$b=NULL
   # plot(hp_new)
    return(hp_new)
  }

  #### normalization of projection vector
  NormalizeProj <- function(projvector){
    p <- nrow(projvector)
    q <- ncol(projvector)
    ss <- sum(projvector[,1]^2)
    normproj <- matrix(stats::rnorm(p*q), ncol=q)
    normproj[,1] <- projvector[,1]/sqrt(ss)

    if(q>1){
      for (j in 1:q){
        temp1 <- 0
        temp2 <- 0
        for (i in 1:p){
          temp1 <- temp1 + projvector[i,j]*normproj[i,j-1]
          temp2 <- temp2 + normproj[i,j-1]*normproj[i,j-1]
        }
        for (i in 1:p){
          normproj[i,j] <- projvector[i,j] - normproj[i,j-1]*temp1/temp2
        }
      }
    }
    return(normproj)
  }

  ############## Projection Pursuit & Find cutoff of clustering ############
  ## find optimal projection vector (p -> q dimension reduction)

  Find.proj<-function(data.now,  PPmethod = "Holes", gamma=0.5,
                      method="kNN",eps=0, 
                      prob=0.2,weight=FALSE,nk=0L,t=1,m=0L,
                      probT=0.05, R=0.1, energy = 0, cooling = 0.999, 
                      TOL = 0.0001, maxiter = 50000L,dataFlag,relCUTOFF=0.1,...){
    

    #data.now <- std.data(data.now)
    tempMean <- colMeans(data.now)
    tempSd <- apply(data.now,2,stats::sd)
    selID = which(tempSd==0)
    if(length(selID)>0){
       tempMean[selID]=0
       tempSd[selID]=1
       data.nowN <- data.now
       data.nowN[,-selID] <- apply(data.now[,-selID],2,function(x) (x-mean(x))/stats::sd(x))
    } else{
       data.nowN <- apply(data.now,2,function(x) (x-mean(x))/stats::sd(x))
    }
#    print(tempSd)
#    print(head(data.nowN))

    data.nowN = as.matrix(data.nowN)
#    print(class(data.nowN))    
    ppopt <- PPclustopt1D(data.nowN, PPmethod = PPmethod, gamma=gamma,
                        method=method,eps=eps, 
                        prob=prob,weight=weight,nk=nk,t=t,m=m,
                        probT=probT, R=R, 
                          energy=energy, cooling=cooling,
                          TOL=TOL, maxiter=maxiter, ...)
    if(ppopt$projbest[1]<0)
      ppopt$projbest = ppopt$projbest*(-1)
    proj.data <- data.nowN%*%ppopt$projbest
    index.best <- ppopt$indexbest


    # P2 : histogram of proj.data
    hist.data <- data.frame(proj.data)
    P2 <- ggplot(data = hist.data, aes(x = proj.data)) +
      geom_histogram( position = "stack") + theme_light()

    best.b <- FindB(proj.data,cutB)
    # add cutoff b to histogram of proj.data
    P2 <- P2 + geom_vline(xintercept = best.b$b)
  #  print(P2)
    L.Flag <- R.Flag <- dataFlag
    L.data <- R.data <- data.now
    if(abs(best.b$rel.dep)<relCUTOFF){
      best.b$b = NULL
      #print(best.b$b)
    }     
    if(!is.null(best.b$b)){
      L.Flag[as.logical(L.Flag)] = proj.data[,1] < best.b$b
      L.data <- L.data[proj.data < best.b$b,]
      R.Flag[as.logical(R.Flag)] <- proj.data[,1] >= best.b$b
      R.data <- R.data[proj.data >= best.b$b,]
    }

    ############ separate data ############

    return(list(Index=index.best,
                proj.vec=ppopt$projbest,
                C=best.b$b,
                L.data=L.data,
                R.data = R.data,
                L.Flag = L.Flag,
                R.Flag = R.Flag,
                P2=P2,
                std.Mean = tempMean,
                std.Sd = tempSd))
  }

  ################### make Cluster Tree ###################
  Tree.construct<-function(data.now,
                           Tree.Struct,
                           id,
                           depth,
                           rep,
                           rep1,
                           rep2,
                           projbest.node,
                           splitCutoff.node,
                           G,
                           histplot.node,
                           final.cluster,
                           cluster.num,
                           split.Mean,
                           split.Sd,
                           dataFlag,relCUTOFF,...) {
    #print(c(id,nrow(data.now),depth))
    ## initialize Tree.Struct
    if(length(Tree.Struct)==0) {
      Tree.Struct<-matrix(1:(2^(max_depth+1)),ncol=1)
      Tree.Struct<-cbind(Tree.Struct,0,0,0,0)
    }

    ## end rule
    if(ifelse(!is.null(max_depth),depth>=max_depth,FALSE)){
      G<-G+1 ####
      Tree.Struct[id,3]<-G ####
      cluster.num <- cluster.num + 1
      #print("FINAL 1")
      #print(depth)
      final.cluster[as.logical(dataFlag)] <- G
      #print(G)
     # print(dataFlag)
     #print(final.cluster)
      return(list(Tree.Struct=Tree.Struct,
                  projbest.node=projbest.node,
                  depth=depth,
                  splitCutoff.node=splitCutoff.node,
                  rep1=rep1,
                  rep2=rep2,
                  G=G,
                  histplot.node=histplot.node,
                  final.cluster=final.cluster,
                  cluster.num=cluster.num,
                  split.Mean = split.Mean,
                  split.Sd= split.Sd,dataFlag = dataFlag,...))
    } else if (nrow(data.now)<minNodeN){
      #print("FINAL 2")
      if (notenough.plot) {gridExtra::grid.arrange(a$P2)}
      G<-G+1 ####
      Tree.Struct[id,3]<-G ####
      cluster.num <- cluster.num + 1
      final.cluster[as.logical(dataFlag)] <- G
      #print(G)
      #print(dataFlag)
      #print(final.cluster)  
      return(list(Tree.Struct=Tree.Struct,
                  projbest.node=projbest.node,
                  depth=depth,
                  splitCutoff.node=splitCutoff.node,
                  rep1=rep1,
                  rep2=rep2,
                  G=G,
                  histplot.node=histplot.node,
                  final.cluster=final.cluster,
                  cluster.num=cluster.num,
                  split.Mean = split.Mean,
                  split.Sd= split.Sd,dataFlag = dataFlag,...))
    } else{
      a <- Find.proj(data.now, PPmethod = PPmethod, gamma=gamma,
                     method=method,eps=eps, 
                     prob=prob,weight=weight,nk=nk,t=t,m=m,
                     probT=probT, R=R, 
                     energy=energy, cooling=cooling,
                     TOL=TOL, maxiter=maxiter,dataFlag=dataFlag,
                     relCUTOFF = relCUTOFF, ...)
     # print(c(print(a$C),nrow(a$L.data),nrow(a$R.data)))
      if (is.null(a$C)|nrow(a$L.data)<minNodeN | nrow(a$R.data)<minNodeN){
       # print("FINAL 3")

        G<-G+1 ####
        Tree.Struct[id,3]<-G ####
        cluster.num <- cluster.num + 1

        final.cluster[as.logical(dataFlag)] <- G
      # print(G)
        #print(dataFlag)
       #print(final.cluster)     
        return(list(Tree.Struct=Tree.Struct,
                    projbest.node=projbest.node,
                    depth=depth,
                    splitCutoff.node=splitCutoff.node,
                    rep1=rep1,
                    rep2=rep2,
                    G=G,
                    histplot.node=histplot.node,
                    final.cluster=final.cluster,
                    cluster.num=cluster.num,
                    split.Mean = split.Mean,
                    split.Sd= split.Sd,dataFlag = dataFlag,...))
      }
      tempMean = a$std.Mean
      tempSd = a$std.Sd


      ## make binary tree
      depth<-depth+1
      Tree.Struct[id,2]<-rep1
      rep1<-rep1+1
      Tree.Struct[id,3]<-rep1
      rep1<-rep1+1
      Tree.Struct[id,4]<-rep2
      rep2<-rep2+1

      splitCutoff.node<-rbind(splitCutoff.node, a$C)
      Tree.Struct[id,5] <- a$Index
      projbest.node <- rbind(projbest.node, t(a$proj.vec))
      histplot.node[[id]] <- a$P2
      split.Mean<- rbind(split.Mean,tempMean)
      split.Sd<- rbind(split.Sd,tempSd)
      b<-Tree.construct(a$L.data, Tree.Struct,
                        id=Tree.Struct[id,2], depth=depth,
                        rep, rep1, rep2, projbest.node, splitCutoff.node,
                        G=G,  histplot.node=histplot.node,
                        final.cluster=final.cluster,
                        cluster.num=cluster.num, 
                        split.Mean = split.Mean,
                        split.Sd= split.Sd,
                        dataFlag=dataFlag*a$L.Flag,relCUTOFF=relCUTOFF,...)
      Tree.Struct<-b$Tree.Struct
      projbest.node<-b$projbest.node
      splitCutoff.node<-b$splitCutoff.node
      rep<-b$rep
      rep1<-b$rep1
      rep2<-b$rep2
      #depth<-b$depth
      G<-b$G ####

      histplot.node<-b$histplot.node
      final.cluster<-b$final.cluster
      cluster.num<-b$cluster.num
      split.Mean <- b$split.Mean
      split.Sd <- b$split.Sd  

      b<-Tree.construct(a$R.data, Tree.Struct,
                        id=Tree.Struct[id,3], depth=depth,
                        rep, rep1, rep2, projbest.node, splitCutoff.node,
                        G=G, histplot.node=histplot.node,
                        final.cluster=final.cluster,
                        cluster.num=cluster.num,
                        split.Mean = split.Mean,
                        split.Sd= split.Sd,
                        dataFlag=dataFlag*a$R.Flag,relCUTOFF=relCUTOFF,...)
      Tree.Struct<-b$Tree.Struct
      projbest.node<-b$projbest.node
      splitCutoff.node<-b$splitCutoff.node
      split.Mean <- b$split.Mean
      split.Sd <- b$split.Sd      
      rep<-b$rep
      rep1<-b$rep1
      rep2<-b$rep2
      depth<-b$depth
      G<-b$G

      histplot.node<-b$histplot.node
      final.cluster<-b$final.cluster
      cluster.num<-b$cluster.num
    }
    return(list(Tree.Struct=Tree.Struct,
                projbest.node=projbest.node,
                splitCutoff.node=splitCutoff.node,
                G=G,
                depth=depth,
                rep=rep,
                rep1=rep1,
                rep2=rep2,
                histplot.node=histplot.node,
                final.cluster=final.cluster,
                cluster.num=cluster.num,
                split.Mean = split.Mean,
                split.Sd= split.Sd,
                dataFlag = dataFlag,...))
  }
  
  dataFlag <- rep(TRUE,nrow(origdata))
  splitCutoff.node<-NULL
  projbest.node<-NULL
  Tree.Struct<-NULL
  final.cluster<-rep(0,nrow(origdata))
  split.Mean <- NULL
  split.Sd <-NULL
  histplot.node<-vector(mode="list", 2^max_depth-1)
  id<-1
  rep1<-2
  rep2<-1
  rep<-1
  cluster.num<-0

  ## make tree
  Tree.final<-Tree.construct(origdata, Tree.Struct, id,
                             depth = 0, rep, rep1, rep2,
                             projbest.node, splitCutoff.node,
                             G=0, histplot.node, final.cluster, cluster.num,
                             split.Mean,split.Sd,
                             dataFlag=dataFlag,relCUTOFF=relCUTOFF,...)
  
  Tree.Struct<-Tree.final$Tree.Struct
  Tree.Struct<-Tree.Struct[-which(Tree.Struct[,3]==0),,drop=FALSE]
  colnames(Tree.Struct)<-c("id","L.node.ID","R.F.node.ID","Coef.ID","Index")
  projbest.node<-Tree.final$projbest.node
  splitCutoff.node<-Tree.final$splitCutoff.node
  histplot.node <- Tree.final$histplot.node
  depth <- Tree.final$depth
  ncluster <- Tree.final$G
  final.cluster <- factor(Tree.final$final.cluster)
  split.Mean <- Tree.final$split.Mean
  split.Sd <- Tree.final$split.Sd
  ## final output
  if(nrow(Tree.Struct)==1){
    Tree.Struct <- NULL
  } 
  treeobj<-list(Tree.Struct=Tree.Struct,
                projbest.node=projbest.node,
                splitCutoff.node=splitCutoff.node,
                origdata=origdata,
                depth=depth,
                histplot.node=histplot.node,
                ncluster=ncluster,
                final.cluster=final.cluster,
                split.Mean=split.Mean,
                split.Sd = split.Sd)

  class(treeobj)<-append(class(treeobj),"PPtreecluster")
  return(treeobj)
}
