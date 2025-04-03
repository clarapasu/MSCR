#include <TMB.hpp>
#include <iostream>
#include <utility>

template<class Type>
Type hazardC(
    Type k0,
    Type k1,
    vector<Type> theta,
    Type t,
    Type z0,
    Type z1,
    Type s0,
    Type s1
){
  Type B = exp(-exp(theta[2]) * t);
  Type mu0 = B * z0 + (1 - B) * s0;
  Type mu1 = B * z1 + (1 - B) * s1;
  Type hz = exp(theta[0]) * exp(-0.5 * ((k0 - mu0) * (k0 - mu0) + (k1 - mu1) * (k1 - mu1)) / (exp(theta[1]) - B * B * exp(theta[1])));
  return hz;
}

template<class Type>
Type halfnormalC(
    Type k0,
    Type k1,
    vector<Type> theta,
    Type s0,
    Type s1){
  Type hn = exp(theta[0]) * exp(-0.5 * ((k0 - s0) * (k0 - s0) + (k1 - s1) * (k1 - s1)) / exp(theta[1]));
  return hn;
}

template<class Type>
Type objective_function<Type>::operator() (){
  PARAMETER_VECTOR(theta);
  
  DATA_MATRIX(trap);
  DATA_MATRIX(df);
  DATA_IVECTOR(dfrows);
  DATA_MATRIX(mesh);
  DATA_SCALAR(endt);
  DATA_INTEGER(n);
  
  int D = mesh.rows();
  int nT = trap.rows();
  
  Type S = 0;
  Type U = 0;
  Type LS = 0;
  Type LL = 0;
  
  vector<Type> theta_vector(theta.size());
  for (int i = 0; i < theta.size(); ++i) {
    theta_vector[i] = theta[i];
  }
  vector<Type> start_hazards(D);
  
  for (int i = 0; i < D; i++) {
    start_hazards[i] = 0;
    for (int j = 0; j < nT; j++) {
      start_hazards[i] += halfnormalC(
        trap(j, 0),     // Keep as Type
        trap(j, 1),     // Keep as Type
        theta_vector,   // Pass the vector<Type> version of theta
        mesh(i, 0),     // Keep as Type
        mesh(i, 1)      // Keep as Type
      );
    }
    U = 1 - exp(-endt * start_hazards[i]);
    S += U;
  }
  
  LS = log(S);
  
  LL = 0;
  int data_startrow = 0;
  int data_endrow = dfrows[0] - 1;
  
  for (int hh = 0; hh < n; hh++) {
    
    matrix<Type> data = df.block(data_startrow, 0, dfrows[hh], 2);
    int m = data.rows();
    vector<Type> t = data.col(0);
    vector<Type> y = data.col(1);
    
    int row_ind = 0;
    // Using a loop to find the first index where y[row_ind] > 0.5
    while (row_ind < m && y[row_ind] <= 0.5) {
      row_ind += 1;
    }
    
    // Setting fcy based on the row index found
    int fcy = (row_ind < m) ? 1 : 0; // Only set fcy if row_ind is valid
    
    Type L = 0;
    for (int j = 0; j < D; j++) {
      Type s0 = mesh(j, 0); // No CppAD::Value needed
      Type s1 = mesh(j, 1); // No CppAD::Value needed
      
      Type capt_time = t[row_ind];
      Type k0 = trap(fcy - 1, 0); // No CppAD::Value needed
      Type k1 = trap(fcy - 1, 1); // No CppAD::Value needed
      
      Type Lj = exp(-(capt_time * start_hazards[j])) * halfnormalC(k0, k1, theta_vector, s0, s1);
      int seen_ind = row_ind;
      
      Type z0 = k0;
      Type z1 = k1;
      
      for (int i = row_ind + 1; i < m; i++) {
        Type total_hazard = 0;
        
        for (int jj = 0; jj < nT; jj++) {
          total_hazard += hazardC(
            trap(jj, 0),      // Keep as Type
            trap(jj, 1),      // Keep as Type
            theta_vector,     // Use the converted vector<Type>
            t[i - 1] - t[seen_ind] + 0.5 * (t[i] - t[i - 1]), 
            z0, z1, 
            s0, s1
          );
        }
        
        Type survival = exp(-(t[i] - t[i - 1]) * total_hazard);
        Lj *= survival;
        
        if (y[i] > 0.5) {
          // Extract the integer value safely using rounding
          int y_index = CppAD::Integer(y[i] + 0.5) - 1; 
          if (y_index < 0) {
            y_index = 0; // Ensure we don't have a negative index
          }
          k0 = trap(y_index, 0);                  // Access trap using the integer index
          k1 = trap(y_index, 1);                  // Access trap using the integer index
          Lj *= hazardC(k0, k1, theta_vector, (t[i] - capt_time), z0, z1, s0, s1);
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
    data_endrow = data_startrow + dfrows[hh + 1] - 1;
  }
  
  LL = LL - n * LS;
  return -LL;
}
