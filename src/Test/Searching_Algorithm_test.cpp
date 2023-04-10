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
TEST(Searching_Algorithm, Three_Moments){

  std::streambuf* orig_buf = std::cout.rdbuf();
  std::cout.rdbuf(NULL);
  
  int num_points = 1000;
  double m = 1.0;

  homedf::Three_Moments::Vector_type V;
  V(0) = 1.3;
  V(1) = 1.25;
  V(2) = 2.65;
  
  auto U     = homedf::Three_Moments::Get_Solution_Vector_from_Primitive(V);
  auto U_ND  = homedf::Three_Moments::Non_Dimensionalize_Solution_Vector(U);
  auto V_ND  = homedf::Three_Moments::Get_Primitive_from_Solution_Vector(U_ND);
  auto F     = homedf::Search_Algorithm<homedf::Three_Moments>(m, V_ND, num_points);  
  auto U_hat = Integrate_Solution_Vector(F, num_points);

  std::cout.rdbuf(orig_buf);
 
  EXPECT_NEAR(U_ND(0), U_hat(0), 1.0e-8);
  EXPECT_NEAR(U_ND(1), U_hat(1), 1.0e-8);
  EXPECT_NEAR(U_ND(2), U_hat(2), 1.0e-8);
}

/////////////////////////////////////////////////////////////////////////
/// Test for Searching Algorithm with Five Moments.
TEST(Searching_Algorithm, Five_Moments){

  std::streambuf* orig_buf = std::cout.rdbuf();
  std::cout.rdbuf(NULL);
  
  int num_points = 1000;
  double m = 1.0;

  homedf::Five_Moments::Vector_type V;
  V(0) = 1.3;
  V(1) = 1.25;
  V(2) = 2.65;
  V(3) = 1.2;
  V(4) = 6.9;
  
  auto U     = homedf::Five_Moments::Get_Solution_Vector_from_Primitive(V);
  auto U_ND  = homedf::Five_Moments::Non_Dimensionalize_Solution_Vector(U);
  auto V_ND  = homedf::Five_Moments::Get_Primitive_from_Solution_Vector(U_ND);
  auto F     = homedf::Search_Algorithm<homedf::Five_Moments>(m, V_ND, num_points);  
  auto U_hat = Integrate_Solution_Vector(F, num_points);

  std::cout.rdbuf(orig_buf);
  
  EXPECT_NEAR(U_ND(0), U_hat(0), 1.0e-8);
  EXPECT_NEAR(U_ND(1), U_hat(1), 1.0e-8);
  EXPECT_NEAR(U_ND(2), U_hat(2), 1.0e-8);
  EXPECT_NEAR(U_ND(3), U_hat(3), 1.0e-8);
  EXPECT_NEAR(U_ND(4), U_hat(4), 1.0e-8);
}

/////////////////////////////////////////////////////////////////////////
/// Test for Searching Algorithm with Seven Moments.
TEST(Searching_Algorithm, Seven_Moments){

  std::streambuf* orig_buf = std::cout.rdbuf();
  std::cout.rdbuf(NULL);
  
  int num_points = 1000;
  double m = 1.0;

  homedf::Seven_Moments::Vector_type V;
  V(0) = 1.2;
  V(1) = 0.15;
  V(2) = 1.1;
  V(3) = 0.2;
  V(4) = 3.5;
  V(5) = 0.1;
  V(6) = 13.0;
  
  auto U     = homedf::Seven_Moments::Get_Solution_Vector_from_Primitive(V);
  auto U_ND  = homedf::Seven_Moments::Non_Dimensionalize_Solution_Vector(U);
  auto V_ND  = homedf::Seven_Moments::Get_Primitive_from_Solution_Vector(U_ND);
  auto F     = homedf::Search_Algorithm<homedf::Seven_Moments>(m, V_ND, num_points);  
  auto U_hat = Integrate_Solution_Vector(F, num_points);

  std::cout.rdbuf(orig_buf);
  
  EXPECT_NEAR(U_ND(0), U_hat(0), 1.0e-8);
  EXPECT_NEAR(U_ND(1), U_hat(1), 1.0e-8);
  EXPECT_NEAR(U_ND(2), U_hat(2), 1.0e-8);
  EXPECT_NEAR(U_ND(3), U_hat(3), 1.0e-8);
  EXPECT_NEAR(U_ND(4), U_hat(4), 1.0e-8);
  EXPECT_NEAR(U_ND(5), U_hat(5), 1.0e-8);
  EXPECT_NEAR(U_ND(6), U_hat(6), 1.0e-8);
}


