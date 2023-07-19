#' @title Visualize PPclusteropt Object
#' @description Visualize the result of projection pursuit optimization
#' @usage PPclustoptViz(PPclustoptOBJ,
#'                      classInfo=NULL)
#' @param PPclustoptOBJ PPclustoptOBJ object. result from PPclustopt1D
#' @param classInfo class information
#' @export
#' @examples
#' data(sampledata)
#' PP.proj.result <- PPclustopt1D(sampledata[,-1])
#' PPclustoptViz(PP.proj.result)
#' PPclustoptViz(PP.proj.result,classInfo=sampledata[,1])
#' @import ggplot2 grid gridExtra PPtreeViz PRIMME

PPclustoptViz<-function(PPclustoptOBJ,classInfo=NULL){
   proj.data<- PPclustoptOBJ$proj.data
   q<-ncol(proj.data)
   p<-ncol(PPclustoptOBJ$origdata)   
   vID <- factor(colnames(PPclustoptOBJ$origdata),
                 levels =colnames(PPclustoptOBJ$origdata) )#1:p
   if(is.null(classInfo)){
     origclass <- rep("",ncol(proj.data))
   } else{
     origclass<-classInfo
   }
   plot.data<-data.frame(proj.data,origclass)
   ..density.. <- NULL
   p1<-ggplot(plot.data,aes(x=proj.data,group=origclass))+
              geom_histogram(aes(y=..density..,fill=origclass))
   coef<-PPclustoptOBJ$projbest#[,1]
   coef.data<-data.frame(vID,coef)
   bin.width<-ifelse(p>100,1,0.1)
   y.max<-max(c(abs(coef.data$coef),1/sqrt(p)))
   y.min<- -y.max
   p2<-ggplot(coef.data,aes(x=vID,y=coef))
   if(p<=10){
      p2<-p2+geom_segment(aes(yend=0,xend=vID,linewidth=1))
   } else{
      p2<-p2+geom_segment(aes(yend=0,xend=vID))
   }       
   p2<-p2+geom_hline(yintercept=0)+
          geom_hline(yintercept=c(-1,1)*1/sqrt(p),
                    col=2,linetype="dashed")+
          ylim(y.min,y.max) +
          xlab("variable ID")+
          ggtitle("Coefficients of Best Projection")+
          theme(legend.position = "none",
                axis.text.x = element_text(angle=45,vjust=1,hjust=1) )
   p3<-ggplot(plot.data,aes(x=proj.data))+
     geom_density()
   gridExtra::grid.arrange(p2,p1,p3,nrow=1)   
}