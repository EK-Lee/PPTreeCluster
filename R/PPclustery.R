#' @title Find Clusters Using PPtreecluster Object
#' @description Find cluster for the test set with the fitted projection pursuit 
#' cluster tree
#' 
#' @usage PPclustery(Tree.result,
#'                   test.data=NULL,
#'                   true.class=NULL,...)
#' @param Tree.result PPtreecluster object 
#' @param test.data  the test dataset
#' @param true.class true class of test dataset if available
#' @param ... arguments to be passed to methods
#' @return cluster predicted cluster
#' @return true.class true class of test dataset from input parameter
#' @return confusion_matrix confusion matrix with cluster and true.class
#' @export
#' @keywords tree clustering
#' @examples
#' data(sampledata)
#' n <- nrow(sampledata)
#' tot <- c(1:n)
#' n.train <- round(n*0.9)
#' train <- sample(tot,n.train)
#' test <- tot[-train]
#' Cluster.result <- PPTreeCluster(sampledata[train,-1])
#' PPclustery(Cluster.result,sampledata[test,-1],sampledata[test,1])

PPclustery<-function(Tree.result,test.data=NULL,true.class=NULL,...) {
  if(is.null(test.data))
    test.data<-Tree.result$origdata
   test.data<-as.matrix(test.data)
   true.class.O <- true.class
   if(!is.null(true.class)){  
      true.class<-as.matrix(true.class); 
      if(nrow(true.class)==1) 
         true.class<-t(true.class)
      if(!is.numeric(true.class)) {
         class.name<-names(table(true.class))
         temp<-rep(0,nrow(true.class))
         for(i in 1:length(class.name))
            temp<-temp+(true.class==class.name[i])*i
         true.class<-temp
      }
   }   

   PP.Classification<-function(Tree.Struct,test.class.index,IOindex,
                               test.class,id,rep){
      if(Tree.Struct[id,4]==0){
         i.class<-test.class
         i.class[i.class>0]<-1
         i.class<-1-i.class
         test.class<-test.class+IOindex*i.class*Tree.Struct[id, 3]
         return(list(test.class=test.class,rep=rep))
      } else{  
         IOindexL<-IOindex*test.class.index[rep,]
         IOindexR<-IOindex*(1-test.class.index[rep,])
         rep<-rep+1
         a<-PP.Classification(Tree.Struct,test.class.index,IOindexL,
                              test.class,Tree.Struct[id,2],rep)
         test.class<-a$test.class
         rep<-a$rep;
         a<-PP.Classification(Tree.Struct,test.class.index,IOindexR,
                              test.class,Tree.Struct[id,3],rep)
         test.class<-a$test.class
         rep<-a$rep
      }
      list(test.class=test.class,rep=rep)
   }
  
   PP.Class.index<-function(class.temp,test.class.index,test.data,
                               Tree.Struct,Alpha.Keep,C.Keep,
                            Split.Mean,Split.Sd,id){
      class.temp<-as.integer(class.temp)
      if(Tree.Struct[id,2]==0){
         return(list(test.class.index=test.class.index,class.temp=class.temp))
      } else{
         t.class<-class.temp 
         t.n<-length(t.class[t.class==0])
         t.index<-sort.list(t.class)
         if(t.n)
            t.index<-sort(t.index[-(1:t.n)])
         t.data<-test.data[t.index,]
         id.proj<-Tree.Struct[id,4]
            
         test.data1 <- as.matrix(t(apply(test.data,1,function(x) 
           (x-Split.Mean[id.proj,])/Split.Sd[id.proj,])))
         
         proj.test<-as.matrix(test.data1)%*%as.matrix(Alpha.Keep[id.proj,])
         proj.test<-as.double(proj.test)
         class.temp<-t(proj.test<C.Keep[id.proj,1]) 
         test.class.index<-rbind(test.class.index,class.temp)
         a<-PP.Class.index(class.temp,test.class.index,test.data,
                           Tree.Struct,Alpha.Keep,C.Keep,
                           Split.Mean,Split.Sd,Tree.Struct[id,2])
         test.class.index<-a$test.class.index
         a<-PP.Class.index(1-class.temp,test.class.index,test.data,
                           Tree.Struct,Alpha.Keep,C.Keep,
                           Split.Mean,Split.Sd,Tree.Struct[id,3])
         test.class.index<-a$test.class.index;
      }
      list(test.class.index=test.class.index,class.temp=class.temp)
   }
    
   if(is.null(Tree.result$Tree.Struct)){
     predict.class <- rep(1,nrow(test.data))
   } else{
   n<-nrow(test.data)
   class.temp<-rep(1,n)
   test.class.index<-NULL
   temp<-PP.Class.index(class.temp,test.class.index,test.data,
                        Tree.result$Tree.Struct,Tree.result$projbest.node,
                        Tree.result$splitCutoff.node,
                        Tree.result$split.Mean, Tree.result$split.Sd,1)
   test.class<-rep(0,n)
   IOindex<-rep(1,n)
   temp<-PP.Classification(Tree.result$Tree.Struct,temp$test.class.index,
                           IOindex,test.class,1,1)
#   if(!is.null(true.class)){
#     temptable <- table(true.class,temp$test.class)
#     if(diff(dim(temptable))<0){
#       true.cluster <- apply(temptable,1,function(x) which.max(x))
#       predict.error<-sum(true.cluster[true.class]!=temp$test.class)
#       predict.error.rate<-mean(true.cluster[true.class]!=temp$test.class)  
#     } else{
#       true.cluster <- apply(temptable,2,function(x) which.max(x))
#       #     predict.error<-sum(true.cluster[true.class]!=temp$test.class)
#       predict.error<-sum(true.cluster[temp$test.class]!=true.class)
#       predict.error.rate<-mean(true.cluster[temp$test.class]!=true.class)        
#     }
# 
#   } else {
#      predict.error<-NA
#      predict.error.rate <- NA
#   }  
#   class.name<-names(table(Tree.result$final.cluster$new))
#   predict.class<-class.name[temp$test.class]
   predict.class <- temp$test.class
   }
   confusion_matrix <- NULL
   if(!is.null(true.class.O)){
     true.class = true.class.O
     confusion_matrix <- table(true.class,predict.class)
   } 

   list(cluster=predict.class,
        true.class = true.class,
        #     predict.error=predict.error, 
        #      predict.error.rate=predict.error.rate,         
        confusion_matrix = confusion_matrix)

}

