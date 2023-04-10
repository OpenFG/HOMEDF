#include<iostream>
#include<iomanip>
#include<string>
#include<sstream>
#include<fstream>
#include<cmath>
#include<math.h>
#include"Eigen/Dense"
#include"../Integration.h"
#include"../Integration_Helper.h"
#include"../Distribution_Function.h"
#include"../Searching_Algorithm.h"
#include"../Output_Messages.h"
#include"../Write_to_VTK.h"
#include"../Five_Moments.h"
#include"../Seven_Moments.h"
#include"omp.h"
#include"mpi.h"


int main(int argc, char* argv[]){
  
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //MPI stuff
  //////////////////////////////////////////////////////////////////////////
  MPI_Init(&argc, &argv);

  int num_proc;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  homedf::Output_Header();
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  int num_points = 1000;
  double m       = 1.0;
  double sigma   = 0.1;
  
  homedf::Distribution_Function<homedf::Five_Moments> F;

  using Vector_type = typename homedf::Distribution_Function<homedf::Five_Moments>::Vector_type;
  
  double rho = 1.0;
  double u   = 0.0;
  double p   = 1.0;
  double q   = 0.0;
  double r   = 3.0;

  double r_max = 6.0;
  
  int num_sigma     = 100;
  double sigma_min  = 0.01;
  double sigma_max  = 0.99;
  double dsigma     = (sigma_max-sigma_min)/static_cast<double>(num_sigma);
  
  int num_q    = 100;
  double q_min = 0.0;
  double q_max = 5.0;
  double dq    = (q_max-q_min)/static_cast<double>(num_q);

  std::vector<double> x((2*num_q-1)*num_sigma);
  std::vector<double> y((2*num_q-1)*num_sigma);
  std::vector<double> z((2*num_q-1)*num_sigma);

  Vector_type V;
  V(0) = rho;
  V(1) = u;
  V(2) = p;
  V(3) = q;
  V(4) = r;
  
  

  /// Loop over q for every value of sigma ///
  for(int i = 0; i < num_sigma; i++){
    sigma = sigma_min + static_cast<double>(i)*dsigma;
    V(3)  = 0.0;
    V(4)  = homedf::Distribution_Function<homedf::Five_Moments>::r_from_sigma(V, sigma);
    
    F = homedf::Search_Algorithm<homedf::Five_Moments>(m, V, num_points);

    q_max = sqrt(2.0*sigma*sigma-sigma*(3.0-r_max));
    dq    = (q_max-q_min)/static_cast<double>(num_q);

    for(int j = 0; j < num_q; j++){
      V(3) = q_min + static_cast<double>(j)*dq;
      V(4) = homedf::Distribution_Function<homedf::Five_Moments>::r_from_sigma(V, sigma);

      F = homedf::Search_Algorithm<homedf::Five_Moments>(m, V, num_points, F);
      
      double s = homedf::Integrate_Velocity_Moment<5>(F, num_points);

      x[(2*num_q-1)*i + (num_q - 1) + j] = V(3);
      y[(2*num_q-1)*i + (num_q - 1) + j] = V(4);
      z[(2*num_q-1)*i + (num_q - 1) + j] = s;
      
    }


    /// To get negative side /////
    V(3)  = 0.0;
    V(4)  = homedf::Distribution_Function<homedf::Five_Moments>::r_from_sigma(V, sigma);
    
    F = homedf::Search_Algorithm<homedf::Five_Moments>(m, V, num_points);
    
    for(int j = 1; j < num_q; j++){
      V(3) = -1.0*(q_min + static_cast<double>(j)*dq);
      V(4) = homedf::Distribution_Function<homedf::Five_Moments>::r_from_sigma(V, sigma);

      F = homedf::Search_Algorithm<homedf::Five_Moments>(m, V, num_points, F);
      
      double s = homedf::Integrate_Velocity_Moment<5>(F, num_points);
      
      x[(2*num_q-1)*i + (num_q - 1) - j] = V(3);
      y[(2*num_q-1)*i + (num_q - 1) - j] = V(4);
      z[(2*num_q-1)*i + (num_q - 1) - j] = s;
  
    }
    
  }


  /// Output to VTK ///
  int num_x = (2*num_q-1);
  int num_y = num_sigma;
  homedf::Write_2D_scatter(x,y,z,num_x,num_y,"Closing_Flux");
  
  MPI_Finalize();
  
  return 0;

}
