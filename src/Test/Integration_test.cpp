#include<iostream>
#include<iomanip>
#include<string>
#include<sstream>
#include<fstream>
#include<cmath>
#include<math.h>
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
/// Test for Integration with 1 gaussian quadrature.
TEST(Integration, 1_Quad){

  double x_min   = 0.0;
  double x_max   = 2.0;
  int num_points = 5000;
  
  auto Integrand = [&] (const double& x){
    return exp(-x*x+2.0*x-1.0);
  };

  double Integral = homedf::gaussian_integrate<1>(Integrand, x_min, x_max, num_points);
  double Integral_EXACT = sqrt(M_PI)*std::erf(1.0);
  
  EXPECT_NEAR(Integral_EXACT, Integral, 1.0e-8);
}

/////////////////////////////////////////////////////////////////////////
/// Test for Integration with 2 gaussian quadrature.
TEST(Integration, 2_Quad){

  double x_min   = 0.0;
  double x_max   = 2.0;
  int num_points = 50;
  
  auto Integrand = [&] (const double& x){
    return exp(-x*x+2.0*x-1.0);
  };

  double Integral = homedf::gaussian_integrate<2>(Integrand, x_min, x_max, num_points);
  double Integral_EXACT = sqrt(M_PI)*std::erf(1.0);
  
  EXPECT_NEAR(Integral_EXACT, Integral, 1.0e-8);
}

/////////////////////////////////////////////////////////////////////////
/// Test for Integration with 3 gaussian quadrature.
TEST(Integration, 3_Quad){

  double x_min   = 0.0;
  double x_max   = 2.0;
  int num_points = 50;
  
  auto Integrand = [&] (const double& x){
    return exp(-x*x+2.0*x-1.0);
  };

  double Integral = homedf::gaussian_integrate<3>(Integrand, x_min, x_max, num_points);
  double Integral_EXACT = sqrt(M_PI)*std::erf(1.0);
  
  EXPECT_NEAR(Integral_EXACT, Integral, 1.0e-8);
}

/////////////////////////////////////////////////////////////////////////
/// Test for Integration with 4 gaussian quadrature.
TEST(Integration, 4_Quad){

  double x_min   = 0.0;
  double x_max   = 2.0;
  int num_points = 50;
  
  auto Integrand = [&] (const double& x){
    return exp(-x*x+2.0*x-1.0);
  };

  double Integral = homedf::gaussian_integrate<4>(Integrand, x_min, x_max, num_points);
  double Integral_EXACT = sqrt(M_PI)*std::erf(1.0);
  
  EXPECT_NEAR(Integral_EXACT, Integral, 1.0e-8);
}

/////////////////////////////////////////////////////////////////////////
/// Test for Integration with 5 gaussian quadrature.
TEST(Integration, 5_Quad){

  double x_min   = 0.0;
  double x_max   = 2.0;
  int num_points = 50;
  
  auto Integrand = [&] (const double& x){
    return exp(-x*x+2.0*x-1.0);
  };

  double Integral = homedf::gaussian_integrate<5>(Integrand, x_min, x_max, num_points);
  double Integral_EXACT = sqrt(M_PI)*std::erf(1.0);
  
  EXPECT_NEAR(Integral_EXACT, Integral, 1.0e-8);
}


