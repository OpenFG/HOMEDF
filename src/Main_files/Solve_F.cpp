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
#include"../Three_Moments.h"
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
  


  //////////////////////////////////////////////////////////////////////////
  /// Solves the Three Moments Distribution Function
  //////////////////////////////////////////////////////////////////////////
  homedf::Three_Moments::Vector_type V3;
  V3(0) = 1.2;
  V3(1) = 0.15;
  V3(2) = 1.1;
  
  auto V3_ND   = homedf::Three_Moments::Non_Dimensionalize_Primitive_Vector(V3);
  auto F3      = homedf::Search_Algorithm<homedf::Three_Moments>(m, V3_ND, num_points);
  
  auto q       = homedf::Integrate_Random_Moment<3>(F3, num_points);
  auto EValue3 = F3.Eigen_Solver(m, num_points);

  homedf::Output_Vector("EigenValues",EValue3);
  homedf::Output_Scalar("q",q);

  Write_Distribution_Function_to_data(F3, num_points, "Distribution_Function_3");
  //////////////////////////////////////////////////////////////////////////

  
  //////////////////////////////////////////////////////////////////////////
  /// Solves the Five Moments Distribution Function
  //////////////////////////////////////////////////////////////////////////
  homedf::Five_Moments::Vector_type V5;
  V5(0) = 1.2;
  V5(1) = 0.15;
  V5(2) = 1.1;
  V5(3) = 0.5;
  V5(4) = 3.5;

  auto V5_ND   = homedf::Five_Moments::Non_Dimensionalize_Primitive_Vector(V5);
  auto F5      = homedf::Search_Algorithm<homedf::Five_Moments>(m, V5_ND, num_points);
  
  auto s       = homedf::Integrate_Random_Moment<5>(F5, num_points);
  auto EValue5 = F5.Eigen_Solver(m, num_points);

  homedf::Output_Vector("EigenValues",EValue5);
  homedf::Output_Scalar("s",s);

  Write_Distribution_Function_to_data(F5, num_points, "Distribution_Function_5");
  //////////////////////////////////////////////////////////////////////////

  
  //////////////////////////////////////////////////////////////////////////
  /// Solves the Seven Moments Distribution Function
  //////////////////////////////////////////////////////////////////////////
  homedf::Seven_Moments::Vector_type V7;
  V7(0) = 1.2;
  V7(1) = 0.15;
  V7(2) = 1.1;
  V7(3) = 0.5;
  V7(4) = 3.5;
  V7(5) = 0.1;
  V7(6) = 15.0;
  
  auto V7_ND   = homedf::Seven_Moments::Non_Dimensionalize_Primitive_Vector(V7);
  auto F7      = homedf::Search_Algorithm<homedf::Seven_Moments>(m, V7_ND, num_points);

  auto w       = homedf::Integrate_Random_Moment<7>(F7, num_points);
  auto EValue7 = F7.Eigen_Solver(m, num_points);

  homedf::Output_Vector("EigenValues",EValue7);
  homedf::Output_Scalar("w",w);

  Write_Distribution_Function_to_data(F7, num_points, "Distribution_Function_7");
  //////////////////////////////////////////////////////////////////////////


  
  MPI_Finalize();
  
  return 0;

}
