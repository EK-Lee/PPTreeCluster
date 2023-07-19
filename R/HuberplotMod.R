#' @title Huber plot
#' @description  Huber plot for 2D data. Plot index values according to 360 
#' degree angle, and histogram of optimally projected data. Modified for the 
#' projection pursuit index without group information. 
#' @usage HuberplotMod(origdata, 
#'                     group=NULL, 
#'                     PPmethod = "Holes", 
#'                     opt.proj = TRUE, 
#'                     UserDefFtn = NULL,
#'                     plotFlag = TRUE,
#'                     main=NULL,...) 
#' @param origdata 2-dimensional numerical data for Huber plot
#' @param group class information vector of data
#' @param PPmethod method for projection pursuit; "Holes", "Skew", "Comb", 
#'                 "NaturalHermite", "Lpp", "P","UserDefProjData", and "UserDefData"
#' @param opt.proj flag to show the best projection in the plot
#' @param UserDefFtn User defined index function when PPmethod="UserDefProjData" or 
#'                   "UserDefData"
#' @param plotFlag flag to show huberplot
#' @param main Main title of the plot
#' @param ... Arguments to be passed to methods.
#' @return huber plot
#' @references ...
#' @import ggplot2 
#' @export
#' @examples
#'  \donttest{
#'  data(sampledata)
#' HuberplotMod(origdata = sampledata[,c(2,3)],
#'              group = sampledata$species,
#'              PPmethod="Holes",
#'              main=paste("Holes Index"))  
#'              
#' HuberplotMod(origdata = sampledata[,2:3],
#'              PPmethod="UserDefData",
#'              UserDefFtn=HolesIndex1D,
#'              main=paste("Holes Index"))
#' }
#'
#'
HuberplotMod <- function (origdata, 
                          group=NULL, 
                          PPmethod = "Holes", 
                          opt.proj = TRUE, 
                          UserDefFtn = NULL,
                          plotFlag = TRUE,
                          main=NULL,
                          ...) {
  index <- NULL
  best.proj <- NULL
  best.index <- 0
  origdata <- as.matrix(origdata)
  origdata <-apply(origdata,2,function(x) (x-mean(x))/stats::sd(x))
  for (i in 0:360) {
    theta <- pi/180 * i
    proj.data <- matrix(cos(theta) * origdata[, 1] + sin(theta) *
                          origdata[, 2])
    proj <- matrix(c(cos(theta), sin(theta)), ncol = 1)
    if (PPmethod == "Holes") {
      newindex <- HolesIndex1D(origdata, proj = proj,...)
    }  else if (PPmethod == "Skew") {
      newindex <- SkewIndex1D(origdata, proj = proj,...)
    }  else if (PPmethod == "Comb") {
      newindex <- CombIndex1D(origdata, proj = proj,...)
    }  else if (PPmethod == "Lpp") {
      newindex <- LppIndex1D(origdata,proj = proj,...)
    }  else if (PPmethod == "NaturalHermite") {
      newindex <- NHIndex1D(origdata,proj = proj,...)
    }  else if (PPmethod == "P") {
      newindex <- PIndex1D(origdata,proj = proj,...)
    }  else if (PPmethod == "UserDefProjData") {
      newindex <- UserDefFtn(proj.data,...)
    }  else if (PPmethod == "UserDefData") {
      newindex <- UserDefFtn(origdata,proj,...)
    }
    index <- c(index, newindex)
  }
  sel.index <- which(index[1:180]==max(index[1:180]))
  theta.best.all <- pi/180 * (sel.index - 1)
  theta.best <- theta.best.all[1]
  proj.data.best <- matrix(cos(theta.best) * origdata[, 1] +
                             sin(theta.best) * origdata[, 2])
  index.best <- max(index)
  range <- round(max(index) - min(index), 5)
  if (range == 0) {
    PPindex <- rep(4, length(index))
  } else {
    PPindex <- (index - min(index))/range * 2 + 3
  }
  data.circle <- NULL
  data.index <- NULL
  for (i in 1:361) {
    theta <- pi/180 * (i - 1)
    data.index <- rbind(data.index, c(PPindex[i] * cos(theta),
                                      PPindex[i] * sin(theta)))
    data.circle <- rbind(data.circle, c(4 * cos(theta), 4 *
                                          sin(theta)))
  }
  maxdiff <- max(c(diff(range(origdata[, 1])), diff(range(origdata[,2]))))
  orig.scaled <- apply(origdata, 2, function(x) (x - mean(x))/maxdiff *3.5)
  data.cX <- data.circle[, 1]
  data.cY <- data.circle[, 2]
  data.X <- data.index[, 1]
  data.Y <- data.index[, 2]
  plot.data <- data.frame(data.cX, data.cY, data.X, data.Y)
  x <- orig.scaled[, 1]
  y <- orig.scaled[, 2]
 # group <- ifelse(is.null(group),"1",group)
  if(is.null(group))
    group <- rep("1",length(x))
  point.data <- data.frame(x, y, group)
  min.X <- min(unlist(plot.data))
  max.X <- max(unlist(plot.data))
  P1 <- ggplot(data = plot.data, aes(x = data.X, y = data.Y)) +
    geom_path() +
    geom_path(aes(x = data.cX, y = data.cY),linetype = "dashed") +
    geom_point(data=point.data,aes(x=x,y=y,color=group,shape = group)) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) + xlab("") + ylab("") +
    coord_fixed() + theme_bw() 

  if (opt.proj) {
    P1 <- P1 + geom_abline(intercept = 0, 
                           slope = sin(theta.best)/cos(theta.best))
    if (length(theta.best.all) > 1)
      for (i in 2:length(theta.best.all))
        P1 <- P1 +
          geom_abline(intercept = 0,
                      slope = sin(theta.best.all[i])/cos(theta.best.all[i]),
                      linetype = "dashed")
  }
  angle<-NA
  P3 <-ggplot(data = data.frame(angle=1:361,PPindex=index),
              aes(angle,PPindex))+
    geom_line()#+ylim(0,1)
  best.proj.data <- proj.data.best
  group <- group
  hist.data <- data.frame(best.proj.data, group)
  P2 <- ggplot(data = hist.data, aes(x = best.proj.data, group = group)) +
    geom_histogram(aes(fill = group), position = "stack")
  if(sum(group!="1")==0){
    P1 <- P1+theme(legend.position = "none",panel.border = element_blank())
    P2 <- P2+theme(legend.position = "none")
  } else{
    P1 = P1+ theme(panel.border = element_blank())
  }

  if(plotFlag){
    lay <- rbind(c(1,2),c(3,2))
    gridExtra::grid.arrange(P1,P2,P3, layout_matrix=lay,
                            top=textGrob(main))
  } else{
    return(index)
  }
}
