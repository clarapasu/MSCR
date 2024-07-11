// Similar functions as in Fit_Func.R but coded in C++ to increase the speed of the computations
#include <Rcpp.h>
using namespace Rcpp;


double hazardC (
    double k0,
    double k1,
    NumericVector theta,       
    double t,     
    double z0,
    double z1,
    double s0,
    double s1
)
{
  
  double B = exp(-exp(theta[2]) * t);
  double mu0 = B * z0 + (1-B) * s0;
  double mu1 = B * z1 + (1-B) * s1;
  double hz = exp(theta[0]) * exp(-0.5 * ((k0 - mu0) * (k0 - mu0)  + (k1 - mu1) * (k1 - mu1))  / (exp(theta[1]) - B * B * exp(theta[1])));
  return hz;
}

double halfnormalC (
    double k0,
    double k1,
    NumericVector theta,        
    double s0,
    double s1
)
{
  
  double hn = exp(theta[0]) * exp(-0.5 * ((k0 - s0) * (k0 - s0)  + (k1 - s1) * (k1 - s1))  / exp(theta[1]));
  return hn;
}

// [[Rcpp::export]]
double LikelihoodC(
    NumericVector theta,
    NumericMatrix trap,
    NumericMatrix df,
    IntegerVector dfrows,
    NumericMatrix mesh,
    double endt) {
  
  int n = unique(df(_, 2)).length();
  
  int D = mesh.nrow();
  
  int nT = trap.nrow();

  
  double S;
  double U;
  double LS;
  double LL;
  
  NumericVector start_hazards(D);
  
  S = 0;
  for (int i = 0; i < D; i++) {
    start_hazards[i] = 0;
    for (int j = 0; j < nT; j++) {
      start_hazards[i] += halfnormalC(trap(j, 0), trap(j, 1), theta, mesh(i, 0), mesh(i, 1));
    }
    U = 1 - exp(-endt * start_hazards[i]);
    S += U;
  }
  LS = log(S);
  
  LL = 0;
  
  int data_startrow;
  int data_endrow;
  int m;
  NumericMatrix data;
  NumericVector t;
  NumericVector y;
  double capt_time;
  
  data_startrow = 0;
  data_endrow = dfrows[0] - 1;
  
  for (int hh = 0; hh < n; hh++) {
    
    data = df(Range(data_startrow, data_endrow) , Range(0, 1));
    m = data.nrow();
    t = data(_, 0);
    y = data(_, 1);
    
    int seen_ind;
    double z0;
    double z1;
    double k0;
    double k1;
    double s0;
    double s1;
    double total_hazard;
    double survival;
    double L;
    double Lj;
    
    int row_ind = 0;
    int fcy = y[row_ind];
    while (fcy == 0) {
      row_ind += 1;
      fcy = y[row_ind];
    }

    L = 0;
    for (int j = 0; j < D; j++) {
      
      s0 = mesh(j, 0);
      s1 = mesh(j, 1);
      
      capt_time = t[row_ind]; 
      
      k0 = trap(fcy - 1, 0);
      k1 = trap(fcy - 1, 1); 
      
      Lj = exp(-(capt_time * start_hazards[j])) * halfnormalC(k0, k1, theta, s0, s1);
      seen_ind = row_ind;
      
      z0 = k0;
      z1 = k1;
      
      for (int i = (row_ind + 1); i < m; i++) {
  
        total_hazard = 0;
        
        for (int jj = 0; jj < nT; jj++) {
          total_hazard += hazardC(trap(jj, 0), trap(jj, 1), theta, t[i-1] - t[seen_ind] + 0.5 * (t[i] - t[i-1]), z0, z1, s0, s1);
        }
        survival = exp(-(t[i]-t[i-1]) * total_hazard);
        
        Lj *= survival;
        
        if (y[i] > 0.5){
          k0 = trap(y[i]-1, 0); 
          k1 = trap(y[i]-1, 1); 
          Lj *= hazardC(k0, k1, theta, (t[i] - capt_time), z0, z1, s0, s1);
          z0 = k0;
          z1 = k1;
          capt_time = t[i];
          seen_ind = i;
        }
      }
      L += Lj; 
      
    }
    
    LL += log(L);
    
    data_startrow += dfrows[hh];
    data_endrow += dfrows[hh+1];

  }
  
  LL = LL - n * LS;
  return -LL;
}


double hazardCnoMem (
    double k0,
    double k1,
    NumericVector theta,       
    double t,     
    double z0,
    double z1,
    double s0,
    double s1
)
{
  
  double B = exp(-exp(100));
  double mu0 = B * z0 + (1-B) * s0;
  double mu1 = B * z1 + (1-B) * s1;
  double hz = exp(theta[0]) * exp(-0.5 * ((k0 - mu0) * (k0 - mu0)  + (k1 - mu1) * (k1 - mu1))  / (exp(theta[1]) - B * B * exp(theta[1])));
  return hz;
}


// [[Rcpp::export]]
double LikelihoodCnoMem(
    NumericVector theta,
    NumericMatrix trap,
    NumericMatrix df,
    IntegerVector dfrows,
    NumericMatrix mesh,
    double endt) {
  
  int n = unique(df(_, 2)).length();
  
  double A = (max(mesh(_,1)) - min(mesh(_,1))) * (max(mesh(_,0)) - min(mesh(_,0)));
  double a = A / mesh.nrow();
  int D = mesh.nrow();
  
  int nT = trap.nrow();
  
  double S;
  double U;
  double LS;
  double LL;
  
  NumericVector start_hazards(D);
  
  S = 0;
  for (int i = 0; i < D; i++) {
    start_hazards[i] = 0;
    for (int j = 0; j < nT; j++) {
      start_hazards[i] += halfnormalC(trap(j, 0), trap(j, 1), theta, mesh(i, 0), mesh(i, 1));
    }
    U = 1 - exp(-endt * start_hazards[i]);
    S += U;
  }
  LS = log(S * a/A);
  
  LL = 0;
  
  int data_startrow;
  int data_endrow;
  int m;
  NumericMatrix data;
  NumericVector t;
  NumericVector y;
  double capt_time;
  
  data_startrow = 0;
  data_endrow = dfrows[0] - 1;
  
  for (int hh = 0; hh < n; hh++) {
    
    data = df(Range(data_startrow, data_endrow) , Range(0, 1));
    m = data.nrow();
    t = data(_, 0);
    y = data(_, 1);

    int seen_ind;
    double z0;
    double z1;
    double k0;
    double k1;
    double s0;
    double s1;
    double total_hazard;
    double survival;
    double L;
    double Lj;
    
    int row_ind = 0;
    int fcy = y[row_ind];
    while (fcy == 0) {
      row_ind += 1;
      fcy = y[row_ind];
    }
    
    L = 0;
    for (int j = 0; j < D; j++) {
      
      s0 = mesh(j, 0);
      s1 = mesh(j, 1);
      
      capt_time = t[row_ind]; 
      
      k0 = trap(fcy - 1, 0); 
      k1 = trap(fcy - 1, 1); 
      
      Lj = exp(-(capt_time * start_hazards[j])) * halfnormalC(k0, k1, theta, s0, s1);
      seen_ind = row_ind;
      
      z0 = k0;
      z1 = k1;
      
      for (int i = (row_ind + 1); i < m; i++) {
        
        total_hazard = 0;
        
        for (int jj = 0; jj < nT; jj++) {
          total_hazard += hazardCnoMem(trap(jj, 0), trap(jj, 1), theta, t[i-1] - t[seen_ind] + 0.5 * (t[i] - t[i-1]), z0, z1, s0, s1);
        }
        survival = exp(-(t[i]-t[i-1]) * total_hazard);
        
        Lj *= survival;
        
        if (y[i] > 0.5){
          k0 = trap(y[i]-1, 0); 
          k1 = trap(y[i]-1, 1); 
          Lj *= hazardCnoMem(k0, k1, theta, (t[i] - capt_time), z0, z1, s0, s1);
          z0 = k0;
          z1 = k1;
          capt_time = t[i];
          seen_ind = i;
        }
      }
      L += Lj; 
      
    }
    
    LL += log(L * a/A);
    
    data_startrow += dfrows[hh];
    data_endrow += dfrows[hh+1];
    
  }
  
  LL = LL - n * LS;
  return -LL;
}
