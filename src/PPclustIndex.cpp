#include "RcppArmadillo.h"
#include "Rcpp.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]


 // [[Rcpp::export]]
 NumericVector matrix2vector(NumericMatrix m, const bool byrow=false){
 if (byrow){
 Rcout << "warning: default by column\n";
 m = transpose(m);
 }
 return NumericVector(m);
 }
 
 
 NumericVector NormalizeProjV(NumericVector proj){
   int p=proj.length();
   NumericVector normproj(p);
   double ss=0;
   for(int i=0;i<p;i++){
     ss+=proj(i)*proj(i);
   }
   for(int i=0;i<p;i++){
     normproj(i)=proj(i)/sqrt(ss);
   }    
   return normproj; 
 }
 
 
 
 Rcpp::IntegerVector Qfloor(Rcpp::NumericVector x) {
 size_t N=x.size();
 Rcpp::IntegerVector y(N);
 for(size_t i=0; i<N; ++i) y[i]=std::floor(x[i]);
 return y;
 }
 
 Rcpp::IntegerVector Qceiling(Rcpp::NumericVector x) {
 size_t N=x.size();
 Rcpp::IntegerVector y(N);
 for(size_t i=0; i<N; ++i) y[i]=std::ceil(x[i]);
 return y;
 }
 
 Rcpp::NumericVector Qsort(Rcpp::NumericVector x) {
 Rcpp::NumericVector y = Rcpp::clone(x);
 std::sort(y.begin(), y.end());
 return y;
 }
 
 // [[Rcpp::export]]
 Rcpp::NumericVector Quantile(Rcpp::NumericVector x,
 Rcpp::NumericVector probs) {
 // implementation of type 7
 const int n=x.size(), np=probs.size();
 if(n==0) return x;
 if(np==0) return probs;
 /*
  int _np=_probs.size(), np;
  np = _np==0 ? 5 : _np;
  Rcpp::NumericVector probs(np);
  if(_np==0) {
  probs[0]=0.00; probs[1]=0.25; probs[2]=0.50; probs[3]=0.75; probs[5]=1.00;
  }
  Rcpp::IntegerVector lo(Rcpp::sapply(index, Qfloor<double>() )), 
  hi(Rcpp::sapply(index, Qceiling<double>() )); 
  */
 Rcpp::NumericVector index=(n-1.)*probs, y=Qsort(x), x_hi(np), qs(np);
 Rcpp::IntegerVector lo(Qfloor(index)), hi(Qceiling(index));
 
 for(size_t i=0; i<np; ++i) {
   qs[i]=y[lo[i]];
   x_hi[i]=y[hi[i]];
   if((index[i]>lo[i]) && (x_hi[i] != qs[i])) {
     double h;
     h=index[i]-lo[i];
     qs[i]=(1.-h)*qs[i]+h*x_hi[i];
   }
 }
 
 return qs;
 //return Rcpp::wrap(qs);
 }
 
// [[Rcpp::export(name=".HolesIndex1D")]]

double HolesIndex1D(NumericMatrix origdata,
                    NumericVector proj=NumericVector(0)){
  
  int n=origdata.nrow(),p=origdata.ncol(),p1=proj.size();
  double index=0;
  double pi = 3.141592;
  NumericVector projdata(n);
  if(p1!=p||p1==1){
    projdata=origdata(_,0);
  } else{
    for(int i=0;i<n;i++){
      for(int k=0;k<p;k++){
        projdata(i)+=origdata(i,k)*proj(k);
      }
    }
  }
  double sproj = 0;
  double projMean = mean(projdata);
  double projSD = sd(projdata);  
  for(int i=0;i<n;i++){
    projdata(i) = (projdata(i)-projMean)/projSD;     
    sproj = sproj+exp(-0.5*projdata(i)*projdata(i))/sqrt(2*pi);
  }
 /* index = ((1/sqrt(2*pi)-sproj/n)*sqrt(2*pi)-0.2)*100/5;*/
 index = (1/sqrt(2*pi)-sproj/n)*sqrt(2*pi)/(1-exp(-0.5));
  return index;
}


// [[Rcpp::export(name=".SkewIndex1D")]]


double SkewIndex1D(NumericMatrix origdata,
                   NumericVector proj=NumericVector(0)){
  
  int n=origdata.nrow(),p=origdata.ncol(),p1=proj.size();
  double index=0;
  double pi = 3.141592;
  double r=0.838;
  double a, b,maxI;
  a=sqrt((1-r)/r);
  b = -sqrt(r/(1-r));
  maxI = r*a*exp(-a*a*0.5)/sqrt(2*pi)-(1-r)*b*exp(-b*b*0.5)/sqrt(2*pi);
  NumericVector projdata(n);
  if(p1!=p||p1==1){
    projdata=origdata(_,0);
  } else{
    for(int i=0;i<n;i++){
      for(int k=0;k<p;k++){
        projdata(i)+=origdata(i,k)*proj(k);
      }
    }
  }
  double projMean = mean(projdata);
  double projSD = sd(projdata);
  double sproj = 0; 
  for(int i=0;i<n;i++){
    projdata(i) = (projdata(i)-projMean)/projSD; 
    sproj = sproj + projdata(i)*exp(-0.5*projdata(i)*projdata(i))/sqrt(2*pi);
  }
  index = abs(sproj)/n;
  return index/maxI;  
}

// [[Rcpp::export(name=".CombIndex1D")]]

double CombIndex1D(NumericMatrix origdata,
                 NumericVector proj=NumericVector(0),
                 double gamma=0.5){
  
  int n=origdata.nrow(),p=origdata.ncol(),p1=proj.size();
  double SKindex=0,Holesindex=0,index=0;
  double pi = 3.141592;
  NumericVector projdata(n);
  if(p1!=p||p1==1){
    projdata=origdata(_,0);
  } else{
    for(int i=0;i<n;i++){
      for(int k=0;k<p;k++){
        projdata(i)+=origdata(i,k)*proj(k);
      }
    }
  }
  double projMean = mean(projdata);
  double projSD = sd(projdata);
  double sprojH = 0,sprojS=0; 
  double ga=0.838;
  double a, b,maxI;
  a=sqrt((1-ga)/ga);
  b = -sqrt(ga/(1-ga));
  maxI = ga*a*exp(-a*a*0.5)/sqrt(2*pi)-(1-ga)*b*exp(-b*b*0.5)/sqrt(2*pi);
  
  for(int i=0;i<n;i++){
    projdata(i) = (projdata(i)-projMean)/projSD; 
    sprojS = sprojS + projdata(i)*exp(-0.5*projdata(i)*projdata(i))/sqrt(2*pi);
    sprojH = sprojH + exp(-0.5*projdata(i)*projdata(i))/sqrt(2*pi);   
   
  }
  Holesindex = (1/sqrt(2*pi)-sprojH/n)*sqrt(2*pi)/(1-exp(-0.5));
  //             ((1/sqrt(2*pi)-sprojH/n)*sqrt(2*pi)-0.2)*100/5;//
  //Holesindex = (1/sqrt(2*pi)-sprojH/n)*sqrt(2*pi)*3;
  //SKindex = abs(sprojS)/n*1000/3;
  SKindex = (abs(sprojS)/n)/maxI;
  //return index/maxI;  
  index = gamma * Holesindex+(1-gamma)*SKindex;
  return index;
}




// [[Rcpp::export(name=".LppIndex1D")]]

double LppIndex1D(NumericMatrix origdata,
                  NumericVector proj=NumericVector(0),
                  std::string  method="kNN",
                  double eps=0,
                  NumericVector probV=NumericVector(0.05),
                  bool weight=false, 
                  int nk=0,double t=1.0){
  
  int n=origdata.nrow(),p=origdata.ncol(),p1=proj.size();  
  if(nk == 0)
    nk = floor(probV(0)*n);

  for(int kk=0;kk<p;kk++){
    double meanO = mean(origdata.column(kk));
    double sdO = sd(origdata.column(kk));     
    origdata.column(kk) = (origdata.column(kk)-meanO)/sdO;
  }

  NumericVector projdata(n);
  for(int i=0;i<n;i++){
      for(int kk=0;kk<p;kk++){
        projdata(i)+=origdata(i,kk)*proj(kk);
      }
  }
  
  double projMean = mean(projdata);
  double projSD = sd(projdata);
  for(int i=0;i<n;i++){
    projdata(i) = (projdata(i)-projMean)/projSD; 
  }
  

  NumericMatrix DISTorig(n,n),DISTproj(n,n);
  
  /* kNN: 1, dist: 2 */
  double result = 0;
  NumericVector selID;
  
  if(method=="kNN"){
    for(int i=0;i<n;i++){
      for(int j=0; j<n; j++){
        if(i!=j){
          DISTorig(i,j)=sum((origdata.row(i)-origdata.row(j))*
            (origdata.row(i)-origdata.row(j)));
          DISTproj(i,j)=
            (projdata(i)-projdata(j))*(projdata(i)-projdata(j));   
        }     
      } 

      NumericVector y = DISTorig.row(i);      
      std::sort(y.begin(),y.end());
      selID = 
        Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(
            arma::sort_index(arma::mat(DISTorig.begin(), 
                                       DISTorig.nrow(), 
                                       DISTorig.ncol(), 
                                       false).row(i)))); 
      for(int jj=1;jj<=nk; jj++){
        if(weight){
          result = result + DISTproj(i,selID(jj))*exp(-DISTorig(i,selID(jj))/t);
        } else{
          result = result + DISTproj(i,selID(jj));
       }
      }
    }
    result = result/nk/n;
  } else if(method=="dist"){
    
    for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
        if(i!=j){
          DISTorig(i,j)=sum((origdata.row(i)-origdata.row(j))*
            (origdata.row(i)-origdata.row(j)));
          DISTproj(i,j)=
            (projdata(i)-projdata(j))*(projdata(i)-projdata(j));   
        }     
      }
    }
    
    if(eps==0){
      NumericVector m = Rcpp::as<Rcpp::NumericVector>(DISTorig);
      eps = Quantile(m,probV)(0); 
      /*# eps = sqrt(eps);*/
    }   
    int count = 0;
    for(int i=0; i<n; i++){
      for(int jj=(i+1); jj<n; jj++){
        if(DISTorig(i,jj)<eps){
          count++;
          if(weight){
            result = result + DISTproj(i,jj)*exp(-DISTorig(i,jj)/t);
          } else{
            result = result + DISTproj(i,jj);
          }
        }
     }
    }
    result = result/count;
  }
  
  double index;
  
  if(result !=0){
    index = (1.0/result);
  } else{
    index = 0.0;
  }

  return(index);
}

double factorial(int num){
   if (num <= 1) return 1;
   return num * factorial(num - 1);
}
 
 // [[Rcpp::export(name=".NHIndex1D")]]
 
double NHIndex1D(NumericMatrix origdata,
                    NumericVector proj=NumericVector(0),
                    int m=0){
   
   int n=origdata.nrow(),p=origdata.ncol(),p1=proj.size();
   double pi = 3.141592;
   
   NumericVector projdata(n);
   if(p1!=p||p1==1){
     projdata=origdata(_,0);
   } else{
     for(int i=0;i<n;i++){
       for(int k=0;k<p;k++){
         projdata(i)+=origdata(i,k)*proj(k);
       }
     }
   }
   double projMean = mean(projdata);
   double projSD = sd(projdata);
   for(int i=0;i<n;i++){
     projdata(i) = (projdata(i)-projMean)/projSD;
   } 
   NumericMatrix M(n,m+1);
   NumericVector x = projdata;
   for(int i=0; i<=m; i++){
     NumericVector temp(n);
     if(i==0){
       for (int kk = 0; kk < n; kk++) 
         temp(kk) = 1;
     } else if(i==1){
        temp =x;
     } else if(i==2){
        temp =x*x  -1;
     } else if(i==3){
        temp =x*x*x  -3*x;
     } else if(i==4){
        temp =x*x*x*x  -6*x*x+ 3;
     } else if(i==5){
        temp =x*x*x*x*x -10*x*x*x+ 15*x;
     } else if(i==6){
        temp =x*x*x*x*x*x -15*x*x*x*x+ 45*x*x-15;
     } else if(i==7){
        temp =x*x*x*x*x*x*x -21*x*x*x*x*x+105*x*x*x -105*x;
     } else if(i==8){
        temp =x*x*x*x*x*x*x*x -28*x*x*x*x*x*x+
         210*x*x*x*x -420*x*x+105;
     } else if(i==9){
        temp =x*x*x*x*x*x*x*x*x -36*x*x*x*x*x*x*x+
         378*x*x*x*x*x-1260*x*x*x+945*x;
     } else if(i==10){
        temp =x*x*x*x*x*x*x*x*x*x-45*x*x*x*x*x*x*x*x+
         630*x*x*x*x*x*x-3150*x*x*x*x+4725*x*x-945;
     } 
     M(_,i) = temp;
   }

   NumericVector phix =  exp(-0.5*projdata*projdata)/sqrt(2*pi);
   double result = 0,a,b,numer,denom;
   for(int ii=0; ii<=m; ii++){
      a = mean(M(_,ii)*phix);
     if(ii%2==1){
        b=0;
     } else{
       int i=ii/2;
       numer =  pow(-1,i)*sqrt(factorial(2*i));
       denom = sqrt(pi)*factorial(i)*(pow(2,(2*i+1)));
       b=numer/denom;
     }
     result = result + (a-b)*(a-b);
   } 
   result = sqrt(result)*(10-m);
   return result;
 }
 
 // [[Rcpp::export(name=".PIndex1D")]]
 
 double PIndex1D(NumericMatrix origdata,
                 NumericVector proj=NumericVector(0),
                 NumericVector probT=NumericVector(0.2), 
                 NumericVector R=NumericVector(0.01)){
   
   int n=origdata.nrow(),p=origdata.ncol(),p1=proj.size();
   double index=0;
   double pi = 3.141592;
   NumericVector projdata(n);
   if(p1!=p||p1==1){
     projdata=origdata(_,0);
   } else{
     for(int i=0;i<n;i++){
       for(int k=0;k<p;k++){
         projdata(i)+=origdata(i,k)*proj(k);
       }
     }
   }

   double LB = Quantile(projdata,probT)(0); 
   double UB = Quantile(projdata,1-probT)(0);  
   double m=0,m2=0;
   int nn=0;
   for(int i=0; i<n; i++){
     if(projdata(i)>=LB && projdata(i)<=UB){
       nn++;
       m += projdata(i);
       m2 += projdata(i)* projdata(i);
     }
   }
   double trimSD = sqrt((m2-m*m/nn)/(nn-1));
   int distN = n*(n-1)/2;
   NumericVector distdata(distN);
   int count=0;
   for(int i=0;i<n;i++){
     for(int j=0; j<i; j++){
         distdata(count)=abs(projdata(i)-projdata(j)); 
       count++;
       }     
     }

   double cutoff = Quantile(distdata, R)(0);
   double mm=0;
   count = 0;
    for(int i=0; i<distN; i++){
       if(distdata(i)<cutoff){
         mm += distdata(i);
         count++;
       }
    }
    double localDEN = cutoff-mm/count;
    index = trimSD*localDEN;
   return index;
 }
 
 // [[Rcpp::export(name=".PPclustopt1D")]]
List PPclustopt1D(NumericMatrix origdata, 
            std::string PPmethod="Holes", double gamma=0.9,
            std::string method="kNN",double eps=0, 
            NumericVector probV=NumericVector(0.2),bool weight=false, 
            int nk=0,double t=1.0, int m=0,
            NumericVector probTV=NumericVector(0), 
            NumericVector RV=NumericVector(0),
            double energy=0,double cooling=0.999, 
            double TOL=0.0001,int maxiter=50000){
   if(nk == 0)
      nk = floor(probV(0)*origdata.nrow());
   if(nk==0)
     nk = 1;
   int n=origdata.nrow(),p=origdata.ncol();
   Environment base("package:base");
   Function table=base["table"];
   GetRNGstate();
   NumericVector projbest(p);
   double indexbest=0,newindex=0;
   if(PPmethod=="Comb"){
     for(int i=0; i<p; i++){
       NumericVector tempproj(p);
       tempproj(i)=1;
       newindex = CombIndex1D(origdata,tempproj,gamma);  
       if(newindex > indexbest) 
         indexbest = newindex;
     } 
   } else if(PPmethod=="Skew"){
     for(int i=0; i<p; i++){
       NumericVector tempproj(p);
       tempproj(i)=1;
       newindex = SkewIndex1D(origdata,tempproj);  
       if(newindex > indexbest) 
         indexbest = newindex;
     }      
   } else if(PPmethod=="Holes"){
     for(int i=0; i<p; i++){
       NumericVector tempproj(p);
       tempproj(i)=1;
       newindex=HolesIndex1D(origdata,tempproj);  
       if(newindex > indexbest) 
         indexbest = newindex;
     } 
   } else if(PPmethod=="Lpp"){
     for(int i=0; i<p; i++){
       NumericVector tempproj(p);
       tempproj(i)=1;
       newindex=LppIndex1D(origdata,tempproj,method,eps,
                            probV,weight, nk,t);   
       if(newindex > indexbest) 
         indexbest = newindex;
     } 
   } else if(PPmethod=="NaturalHermite"){
     for(int i=0; i<p; i++){
       NumericVector tempproj(p);
       tempproj(i)=1;
       newindex=NHIndex1D(origdata,tempproj,m);  
       if(newindex > indexbest) 
         indexbest = newindex;
     } 
   } else if(PPmethod=="P"){
     for(int i=0; i<p; i++){
       NumericVector tempproj(p);
       tempproj(i)=1;
       newindex=PIndex1D(origdata,tempproj,probTV,RV);  
       if(newindex > indexbest) 
         indexbest = newindex;
     } 
   }      
   if(energy==0)
     energy=abs(1-indexbest);
   if(energy>1)
     energy = 1;
     projbest=rnorm(p);

   projbest=NormalizeProjV(projbest);  
  
   if(PPmethod=="Comb"){
     indexbest = CombIndex1D(origdata,projbest,gamma);  
   } else if(PPmethod=="Skew"){
     indexbest = SkewIndex1D(origdata,projbest);    
   } else if(PPmethod=="Holes"){
     indexbest = HolesIndex1D(origdata,projbest);  
   } else if(PPmethod=="Lpp"){
     indexbest = LppIndex1D(origdata,projbest,method,eps,
                            probV,weight, nk,t);     
   } else if(PPmethod=="NaturalHermite"){
     indexbest = NHIndex1D(origdata,projbest,m);
   }  else if(PPmethod=="P"){
     indexbest = PIndex1D(origdata,projbest,probTV,RV);  
   }    
   double temp=1;
   int kk=0;  
   double diff=100;
   while(fabs(diff)>TOL&&kk<maxiter){
     double tempp=energy/log(kk+2.0);
     if(kk>1000) {
       temp=temp*cooling;
     } else {
       temp=temp*cooling*cooling;
     }   
     
     NumericVector projnew(p);
     projnew=temp*rnorm(p)+projbest;

     projnew=NormalizeProjV(projnew);  
     
   
     if(PPmethod=="Comb"){
       newindex = CombIndex1D(origdata,projnew,gamma);  
     } else if(PPmethod=="Skew"){
       newindex = SkewIndex1D(origdata,projnew);    
     } else if(PPmethod=="Holes"){
       newindex = HolesIndex1D(origdata,projnew);  
     } else if(PPmethod=="Lpp"){
       newindex = LppIndex1D(origdata,projnew,method,eps,
                             probV,weight, nk,t);     
     } else if(PPmethod=="NaturalHermite"){
       newindex = NHIndex1D(origdata,projnew,m);
     }  else if(PPmethod=="P"){
       newindex = PIndex1D(origdata,projnew,probTV,RV);  
     }      
     
    NumericVector prob=runif(1);
     double difft=newindex-indexbest;
     double e=exp(difft/tempp);
     if(e>1){
       for(int i=0;i<p;i++){
           projbest(i)=projnew(i);
       }
       indexbest=newindex;    
       diff=difft;
     } else if(prob[0]<e&&difft>energy){
       for(int i=0;i<p;i++){
           projbest(i)=projnew(i);
       }
       indexbest=newindex;    
       diff=difft; 
     }
     kk++;
   }
   PutRNGstate();
   return Rcpp::List::create(Rcpp::Named("indexbest")=indexbest,
                             Rcpp::Named("projbest")=projbest,
                             Rcpp::Named("origdata")=origdata);
}
 