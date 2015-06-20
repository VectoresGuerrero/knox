//Ported from R code at http://www.niph.go.jp/soshiki/gijutsu/download/Rfunctions/func/KnoxM.test
//#include <Rcpp.h>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//http://stackoverflow.com/questions/24253228/performance-of-r-statssd-vs-armastddev-vs-rcpp-implementation
double armaSD2(const arma::colvec & inVec) { return arma::stddev(inVec); }

//' @rdname knox
//' @title Original Knox test
//' @description Original Knox test
//' @param x Longitude, should be projected to a planar system.
//' @param y Latitude, should be projected to a planar system.
//' @param time time column.
//' @param sigma a distance in space.
//' @param tau a distance in time.
//' @return A knox statistic, which is the number of pairs of points found in a given space-time distance.
//' @references \strong{Knox, E. (1964)}. \emph{The detection of space-time interactions.} Journal of the Royal Statistical Society. Series C (13(1), 25-30.
//' 
//' \strong{Tango, T. (2010)}. \emph{Statistical methods for disease clustering. Springer.}
//' @export
// [[Rcpp::export]]
List knox(NumericVector x, NumericVector y, NumericVector time, double sigma, double tau){
  int nx = x.size();
  NumericMatrix sdist(nx,nx);
  NumericMatrix tdist(nx,nx);
  NumericMatrix as(nx,nx);
  NumericMatrix at(nx,nx);
  //  NumericMatrix atmont(nx,nx);
  double s1,obs;
  
  for(int i=0;i<nx;++i){
    for (int j=0; j<nx; j++){
      sdist(i,j)= sqrt(pow(x[j]-x[i],2)+pow(y[j]-y[i],2));
      tdist(i,j)=abs(time[i]-time[j]);
      if (i !=j){
        as(i,j) = (sdist(i,j)<= sigma) ? 1 : 0;
        at(i,j) = (tdist(i,j)<= tau) ? 1 : 0;
      }
    }
  }
  s1 = 0;
  for(int i=0;i<nx;++i){
    for(int j=0;j<nx;++j){
      s1 += as(i,j)*at(i,j);
    }
  }
  obs = s1/2;
  return List::create(
    _["knox"]=obs
    );
}

//' @rdname knox_simulation 
//' @title Knox test with Monte Carlo simulation
//' @description Performing Knox statistics with Monte Carlo simulation
//' @param x Longitude, should be projected to a planar system.
//' @param y Latitude, should be projected to a planar system.
//' @param time time column.
//' @param sigma a distance in space.
//' @param tau a distance in time.
//' @param perm Number of simulations. Default is 999
//' @return \item{knox}{Knox statistic, which is the number of pairs of points found in a given space-time distance.} 
//' \item{p_value}{p-value calculated from simulation} 
//' \item{RR}{Relative Risk - calculated by observed value (Knox statistics) divided by mean of simulated values}
//' @export
// [[Rcpp::export]]
List knox_mc(NumericVector x, NumericVector y, NumericVector time, double sigma, double tau,int perm=999){
  RNGScope scope;//Initiate to get sample function from R
  Environment base("package:base");
  Function sample = base["sample"];
  int nx = x.size();
  NumericMatrix sdist(nx,nx);
  NumericMatrix tdist(nx,nx);
  NumericMatrix as(nx,nx);
  NumericMatrix at(nx,nx);
  int perm1 = perm + 1; //Normally 999 interations + 1 which is real obs Knox value to get 1000. Since then p_value is calculated easier
  double ktmont[perm1]; //Create an C++ array for simulated Knox value
  double s1,obs,p_value;
  double sum_perm=0;
  NumericVector timeR(nx);
  
  for(int i=0;i<nx;++i){
    for (int j=0; j<nx; j++){
      sdist(i,j)= sqrt(pow(x[j]-x[i],2)+pow(y[j]-y[i],2));
      tdist(i,j)=abs(time[i]-time[j]);
      if (i!=j){
      as(i,j) = (sdist(i,j)<= sigma) ? 1 : 0; //this saved my life
      at(i,j) = (tdist(i,j)<= tau) ? 1 : 0;  
      }
    }
  }
  s1 = 0;
  for(int i=0;i<nx;++i){
    for(int j=0;j<nx;++j){
      s1 += as(i,j)*at(i,j);
    }
  }
  obs = s1/2;
  //Monte Carlo Simulation
  for(int k=0;k<perm;++k){
    std::fill(at.begin(),at.end(),0);
    timeR = sample(time);
    for(int i=0;i<nx;++i)
    for (int j=0;j<nx;++j){
      tdist(i,j)=abs(timeR[i]-timeR[j]);
      if (i!=j){
        at(i,j) = (tdist(i,j)<= tau) ? 1 : 0;  
      }
      
    }
    s1=0;
    for(int i=0;i<nx;++i)
    for(int j=0;j<nx;++j){
      s1 +=as(i,j)*at(i,j);
    }
    ktmont[k]=s1/2;
    sum_perm += ktmont[k];
  }//End of Monte Carlo Simulation
  
  ktmont[perm1]=obs;
  
  double countm=0;
  for(int i=0;i<perm;i++){
    if(ktmont[i]>=obs){
      ++countm;
    }
  }
  p_value = countm/perm1;
  double mean_sum_perm = sum_perm/perm;
  double rr = obs/mean_sum_perm;
  return List::create(
    _["knox"]=obs,
    _["p_value"]=p_value,
    _["RR"]=rr
    );
}


// [[Rcpp::export]]
List iknox_mc(NumericVector x, NumericVector y, NumericVector time, double sigma, double tau,int perm=999){
  RNGScope scope;//Initiate to get sample function from R
  Environment base("package:base");
  Function sample = base["sample"];
  int nx = x.size();
  NumericMatrix sdist(nx,nx);
  NumericMatrix tdist(nx,nx);
  NumericMatrix as(nx,nx);
  NumericMatrix at(nx,nx);
  int perm1 = perm + 1;
  NumericVector ktmont(perm1);
  double s1,obs,p_value;
  double sum_perm=0;
  NumericVector timeR(nx);
  
  for(int i=0;i<nx;++i){
    for (int j=0; j<nx; j++){
      sdist(i,j)= sqrt(pow(x[j]-x[i],2)+pow(y[j]-y[i],2));
      tdist(i,j)=abs(time[i]-time[j]);
      if(i!=j){
        as(i,j) = (sdist(i,j)<= sigma) ? 1 : 0;
        at(i,j) = (tdist(i,j) == tau) ? 1 : 0;        
      }
    }
  }
  s1 = 0;
  for(int i=0;i<nx;++i){
    for(int j=0;j<nx;++j){
      s1 += as(i,j)*at(i,j);
    }
  }
  obs = s1/2;
  //Monte Carlo Simulation
  for(int k=0;k<perm;++k){
    std::fill(at.begin(),at.end(),0);
    timeR = sample(time);
    for(int i=0;i<nx;++i)
    for (int j=0; j<nx;++j){
      tdist(i,j)=abs(timeR[i]-timeR[j]);
      if(i!=j){
        at(i,j) = (tdist(i,j)==tau) ? 1 : 0;         
      }
    }
    s1=0;
    for(int i=0;i<nx;++i)
    for(int j=0;j<nx;++j){
      s1 += as(i,j)*at(i,j);
    }
    ktmont(k)=s1/2;
    sum_perm += ktmont(k);
  }//End of Monte Carlo Simulation
  
  ktmont(perm)=obs;
  
  double countm=0;
  for(int i=0;i<perm1;i++){
    if(ktmont(i)>=obs){
      ++countm;
    }
  }
  //P value
  p_value = countm/perm1;
  //calculate sd using arma lib as its results faster than standard one.
  double sd = armaSD2(ktmont);
  double mean_sum_perm = sum_perm/perm;
  //Knox statistics standardized
  double z = (obs-sd)/mean_sum_perm;
  double rr = obs/mean_sum_perm;
  
  return List::create(
    _["IKT"]=obs,
    _["Z_score"]=z,
    _["p_value"]=p_value,
    _["RR"]=rr
    );
}

//' @rdname st_links
//' @title space-time link
//' @description Create space-time links between two points from a given space and time distance.
//' @param x Longitude - should be projected to a planar system.
//' @param y Latitude - should be projected to a plannar system.
//' @param time time
//' @param ds a distance in space.
//' @param dt a distance in time.
//' @return A data frame with \item{Xo, Yo}{ X, Y of the starting point}
//' \item{Xd, Yd}{ X, Y of the ending point}
//' @export
//[[Rcpp::export]]
DataFrame st_link(NumericVector x, NumericVector y, NumericVector time, double ds, double dt){
  int nx = x.size();
  NumericMatrix sdist(nx,nx);
  NumericMatrix tdist(nx,nx);
  std::vector<double> xo,yo,xd,yd;
  //Calculate space and time distance as a matrix
  for(int i=0;i<nx;++i)
  for (int j=0; j<nx;++j){
    sdist(i,j)= sqrt(pow(x[j]-x[i],2)+pow(y[j]-y[i],2));
    tdist(i,j)=abs(time[i]-time[j]);
    if(j>i && sdist(i,j) <= ds && tdist(i,j) <=dt){
      xo.push_back(x(i));
      yo.push_back(y(i));
      xd.push_back(x(j));
      yd.push_back(y(j));
    }
  }
  return DataFrame::create(Named("Xo")=xo,Named("Yo")=yo,
  Named("Xd")=xd,Named("Yd")=yd);
}
