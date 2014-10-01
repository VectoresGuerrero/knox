//Ported from R code at http://www.niph.go.jp/soshiki/gijutsu/download/Rfunctions/func/KnoxM.test
//#include <Rcpp.h>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//http://stackoverflow.com/questions/24253228/performance-of-r-statssd-vs-armastddev-vs-rcpp-implementation
double armaSD2(const arma::colvec & inVec) { return arma::stddev(inVec); }

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
      //      if(sdist(i,j)<= sigma){
      //        as(i,j)=1;
      //      }
      if (i !=j){
        as(i,j) = (sdist(i,j)<= sigma) ? 1 : 0;
        at(i,j) = (tdist(i,j)<= tau) ? 1 : 0;
      }
      
      //      if(tdist(i,j)<= tau){
      //        at(i,j)=1;
      //      }
    }
  }
  //    for(int i=0;i<nx;i++){
  //      for(int j=0;j<nx;j++){
  //        if(sdist(i,j)<= sigma){
  //          as(i,j)=1;
  //        } 
  //        else { 
  //          as(i,j)=0;
  //        }
  //        if(tdist(i,j)<= tau){
  //          at(i,j)=1;
  //        }
  //        else {
  //          at(i,j)=0;
  //        }
  //      }
  //      as(i,i)=0;
  //      at(i,i)=0;
  //    }
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
  //  NumericMatrix atmont(nx,nx);
  int perm1 = perm + 1;
  double ktmont[perm1];
  double s1,obs,p_value;
  double sum_perm=0;
  NumericVector timeR(nx);
  
  for(int i=0;i<nx;++i){
    for (int j=0; j<nx; j++){
      sdist(i,j)= sqrt(pow(x[j]-x[i],2)+pow(y[j]-y[i],2));
      tdist(i,j)=abs(time[i]-time[j]);
//      if(sdist(i,j)<= sigma){
//        as(i,j)=1;
//      }
      if (i!=j){
      as(i,j) = (sdist(i,j)<= sigma) ? 1 : 0; //this saved my life
      at(i,j) = (tdist(i,j)<= tau) ? 1 : 0;  
      }
//      if(tdist(i,j)<= tau){
//        at(i,j)=1;
//      }
    }
  }
  //naive loop and condition
  //  for(int i=0;i<nx;i++){
  //    for(int j=0;j<nx;j++){
  //      if(sdist(i,j)<= sigma){
  //        as(i,j)=1;
  //      } 
  //      else { 
  //        as(i,j)=0;
  //      }
  //      if(tdist(i,j)<= tau){
  //        at(i,j)=1;
  //      }
  //      else {
  //        at(i,j)=0;
  //      }
  //    }
  //    as(i,i)=0;
  //    at(i,i)=0;
  //  }
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
      //      if(tdist(i,j)<= tau){
      //        at(i,j)=1;
      //      }
      if (i!=j){
        at(i,j) = (tdist(i,j)<= tau) ? 1 : 0;  
      }
      
    }
    //    for(int i=0;i<nx;++i){
    //      for(int j=0;j<nx;++j){
    //        if(tdist(i,j)<= tau){
    //        at(i,j)=1;
    //      } else {
    //        at(i,j)=0;
    //      }
    //    }
    //    at(i,i)=0;
    //  }
    s1=0;
    for(int i=0;i<nx;++i)
    for(int j=0;j<nx;++j){
      s1 +=as(i,j)*at(i,j);
    }
    ktmont[k]=s1/2;
    sum_perm += ktmont[k];
  }//End of Monte Carlo Simulation
  
  ktmont[perm]=obs;
  
  double countm=0;
  for(int i=0;i<perm1;i++){
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
    //return ktmont;
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
  //  NumericMatrix atmont(nx,nx);
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
  //calculate sd using arma lib as it results faster than standard one.
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
  //  for (int i=0;i<nx;++i){
  //    for(int j=0;j<nx;++j){
  //      if(j>i && sdist(i,j) <= ds && tdist(i,j) <=dt){
  //        xo.push_back(x(i));        
  //        yo.push_back(y(i));
  //        xd.push_back(x(j));
  //        yd.push_back(y(j));
  //      }
  //    }
  //  }
  return DataFrame::create(Named("Xo")=xo,Named("Yo")=yo,
  Named("Xd")=xd,Named("Yd")=yd);
}
