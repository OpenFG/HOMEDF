#include<iostream>
#include<iomanip>
#include<string>
#include<sstream>
#include<fstream>
#include<cmath>
#include<math.h>
#include<iostream>
#include"Eigen/Dense"
#include"gtest/gtest.h"
#include"../Integration.h"
#include"../Integration_Helper.h"
#include"../Distribution_Function.h"
#include"../Searching_Algorithm.h"
#include"../Output_Messages.h"
#include"../Three_Moments.h"
#include"../Five_Moments.h"
#include"../Seven_Moments.h"


/////////////////////////////////////////////////////////////////////////
/// Test for Searching Algorithm with Three Moments.
TEST(Distribution_Function, Three_Moment_Wavespeeds){

  std::streambuf* orig_buf = std::cout.rdbuf();
  std::cout.rdbuf(NULL);
  
  int num_points = 1000;
  double m = 1.0;

  homedf::Three_Moments::Vector_type V;
  V(0) = 1.3;
  V(1) = 1.25;
  V(2) = 2.65;
  
  auto V_ND  = homedf::Three_Moments::Non_Dimensionalize_Primitive_Vector(V);
  auto F     = homedf::Search_Algorithm<homedf::Three_Moments>(m, V_ND, num_points);  
  auto EV    = (F.Eigen_Solver(m, num_points)).real();

  std::sort(EV.data(), EV.data()+EV.size());
  
  homedf::Three_Moments::Vector_type EV_EXACT;
  EV_EXACT(0) = V_ND(1) - sqrt(3.0*V_ND(2)/V_ND(0));
  EV_EXACT(1) = V_ND(1);
  EV_EXACT(2) = V_ND(1) + sqrt(3.0*V_ND(2)/V_ND(0));
    
  std::cout.rdbuf(orig_buf);
 
  EXPECT_NEAR(EV_EXACT(0), EV(0), 1.0e-8);
  EXPECT_NEAR(EV_EXACT(1), EV(1), 1.0e-8);
  EXPECT_NEAR(EV_EXACT(2), EV(2), 1.0e-8);
}




