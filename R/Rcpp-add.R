#' @title Holes Index Derived From Natural Hermite Index 
#' @description  Find Central Hole distribution
#' @usage HolesIndex1D(origdata,
#'                     proj)
#' @param origdata the original data matrix
#' @param proj projection vector you want to know the index value
#' @return numeric index value
#' @references Cook, D., Buja, A., & Cabrera, J. (1993). Projection pursuit 
#' indexes based on orthonormal function expansions. Journal of Computational 
#' and Graphical Statistics, 2(3), 225-250.
#' @import ggplot2 quantreg optimx
#' @export
#' @examples
#'  \donttest{
#' data(sampledata)   
#' HolesIndex1D(sampledata[,-1],c(1,0,0,0))
#' }
#' @useDynLib PPtreeCluster
#' @importFrom Rcpp evalCpp
#'
HolesIndex1D <- function(origdata, proj) {
  origdata = as.matrix(origdata)
  .HolesIndex1D(origdata, proj)
}


#' @title Skewness Index Derived From Natural Hermite Index 
#' @description  Find skewed distribution
#' @usage SkewIndex1D(origdata,
#'                    proj)
#' @param origdata the original data matrix
#' @param proj projection vector you want to know the index value
#' @return numeric index value
#' @references Cook, D., Buja, A., & Cabrera, J. (1993). Projection pursuit 
#' indexes based on orthonormal function expansions. Journal of Computational 
#' and Graphical Statistics, 2(3), 225-250.
#' @import ggplot2 quantreg optimx
#' @export
#' @examples
#'  \donttest{
#' data(sampledata)
#' SkewIndex1D(sampledata[,-1],c(1,0,0,0))
#' }

SkewIndex1D <- function(origdata, proj) {
  origdata = as.matrix(origdata)
  .SkewIndex1D(origdata, proj)
}


#' @title Combined Index with Holes and Skewness Index
#' @description  Find interesting patterns from data
#' @usage CombIndex1D(origdata,
#'                    proj,
#'                    gamma=0.9)
#' @param origdata the original data matrix
#' @param proj projection vector you want to know the index value
#' @param gamma ratio for combining two indices
#' @return numeric index value
#' @import ggplot2 quantreg optimx
#' @export
#' @examples
#'  \donttest{
#' data(sampledata)
#' CombIndex1D(sampledata[,-1],c(1,0,0,0))
#' }
#' 
CombIndex1D <- function(origdata, proj,gamma=0.9) {
  origdata = as.matrix(origdata)
  .CombIndex1D(origdata, proj,gamma)
}

#' @title Lpp Index (Espezua et al., 2015)
#' @description  Find projections in which close points on the original 
#'               data are also located close after projection
#' @usage LppIndex1D(origdata, 
#'                   proj, 
#'                   method="kNN",
#'                   eps=0,
#'                   prob=0.2,
#'                   weight=FALSE,
#'                   nk=0,
#'                   t=1)
#' @param origdata the original data matrix
#' @param proj projection vector you want to know the index value
#' @param method "kNN" determines neighbors based on count, while "dist" 
#'               determines based on distance.
#' @param eps distance criterion at the "dist" method
#' @param prob ratio for determining the neighborhood for kNN method
#' @param weight If "TRUE", the target function use a weight based on distance; 
#'               if "FALSE", no weight
#' @param nk count criterion at the "kNN" method
#' @param t weight parameter
#' @return numeric index value
#' @references Espezua, S., Villanueva, E., Maciel, C. D., & Carvalho, A. 
#' (2015). A Projection Pursuit framework for supervised dimension reduction of 
#' high dimensional small sample datasets. Neurocomputing, 149, 767-776.
#' @references He, X., & Niyogi, P. (2003). Locality preserving projections. 
#' Advances in neural information processing systems, 16.
#' @import ggplot2 quantreg optimx
#' @export
#' @examples
#'  \donttest{
#'  data(sampledata)
#'  LppIndex1D(sampledata[,-1],c(1,0,0,0))
#' }
#' 
LppIndex1D <- function(origdata, proj, method="kNN",eps=0, 
                       prob=0.2,weight=FALSE,nk=0,t=1){
    #metho=ifelse(method=="kNN",1,2)
    #eps = ifelse(is.null(eps),0,eps)
    #weight = as.numeric(weight)
    nk = ifelse(nk==0,floor(prob*nrow(origdata)),nk)
    origdata = as.matrix(origdata)
  .LppIndex1D(origdata, proj,method,eps,prob,weight,nk,t)
}


#' @title Natural Hermite Index 
#' @description  Fine projections away from normality
#' @usage NHIndex1D(origdata,
#'                  proj,
#'                  m=0)
#' @param origdata the original data matrix
#' @param proj projection vector you want to know the index value
#' @param m order of the index. low order captures the overall structure 
#'          of data, high order finds detailed patterns of data
#' @return numeric index value
#' @references Cook, D., Buja, A., & Cabrera, J. (1993). Projection pursuit 
#' indexes based on orthonormal function expansions. Journal of Computational 
#' and Graphical Statistics, 2(3), 225-250.
#' @import ggplot2 quantreg optimx
#' @export
#' @examples
#'  \donttest{
#'  data(sampledata)
#'  NHIndex1D(sampledata[,-1],c(1,0,0,0))
#' }
#' 
NHIndex1D <- function(origdata, proj,m=0) {
  origdata = as.matrix(origdata)
  .NHIndex1D(origdata, proj,m)
}

#' @title P Index (Friedman et al., 1974)
#' @description  Find a direction that data points are gathered inside the 
#'               clusters, at the same time, separating the clusters
#' @usage PIndex1D(origdata,
#'                 proj,
#'                 probT=0.05,
#'                 R=0.1)
#' @param origdata the original data matrix
#' @param proj projection vector you want to know the index value
#' @param probT the ratio to trim data
#' @param R criterion for determining neighbors
#' @return numeric index value
#' @references Friedman, J. H., & Tukey, J. W. (1974). A projection pursuit 
#' algorithm for exploratory data analysis. IEEE Transactions on computers, 
#' 100(9), 881-890.
#' @import ggplot2 quantreg optimx PPCI
#' @export
#' @examples
#'  \donttest{
#' data(sampledata)
#' PIndex1D(sampledata[,-1],c(1,0,0,0))
#' }
#'
PIndex1D <- function(origdata, proj,probT=0.05, R=0.1) {
  origdata = as.matrix(origdata)
  .PIndex1D(origdata, proj,probT,R)
}

#' @title Optimization for Projection Pursuit Indices 
#' @description
#' Find the 1-dim optimal projection using various projection pursuit indices 
#' with class information
#' @usage PPclustopt1D(origdata,  
#'                     PPmethod = "Holes", 
#'                     gamma=0.9, 
#'                     method="kNN",
#'                     eps=0, 
#'                     prob=0.2,
#'                     weight=FALSE,
#'                     nk=0L,
#'                     t=1,
#'                     m=0L,
#'                     probT=0.05, 
#'                     R=0.1, 
#'                     energy = 0, 
#'                     cooling = 0.999,
#'                     TOL = 0.0001, 
#'                     maxiter = 50000L) 
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
#' @param energy energy parameter
#' @param cooling cooling parameter
#' @param TOL tolerance
#' @param maxiter number of maximum iteration
#' @return indexbest maximum projection pursuit index value
#' @return projbest optimal 1-dim projection matrix
#' @return origdata  original data matrix 
#' @return proj.data projected data onto the best projection 
#' @export
#' @keywords projection pursuit
#' @examples
#' data(sampledata)
#' PP.proj.result <- PPclustopt1D(sampledata[,-1])
#' hist(PP.proj.result$proj.data)
PPclustopt1D <- function(origdata,  PPmethod = "Holes", gamma=0.9,
                       method="kNN",eps=0, 
                       prob=0.2,weight=FALSE,nk=0L,t=1,m=0L,
                       probT=0.05, R=0.1, energy = 0, cooling = 0.999, 
                       TOL = 0.0001, maxiter = 50000L) {
#  stdData = apply(origdata,2,function(x)  (x-mean(x))/sd(x))
  origdata = as.matrix(origdata)
  result<-.PPclustopt1D(origdata, PPmethod, 
                      gamma, method,eps,  prob, weight, nk,t,m,
                      probT,R, energy, cooling, TOL, maxiter)
#  result$origdata <- origdata
  result$proj.data = origdata%*%as.matrix(result$projbest)
  class(result)<-append(class(result),"PPclustoptOBJ")
  return(result)
}
